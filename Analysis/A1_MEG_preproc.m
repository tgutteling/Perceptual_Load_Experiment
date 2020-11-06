function [] = A1_MEG_preproc(sub)
%
% Initial preprocessing of MEGIN MEG data
%
%   Usage:
%     [] = A_preproc(sub)
%
%     Inputs:
%
%     sub - Subject identifier code (string) e.g. 'S01'
%

%% Set up folders

%set path
run('/rds/projects/2017/jenseno-01/Tjerk/set_path');

%add scripts location
addpath('/rds/projects/2017/jenseno-01/Tjerk/Load2/scripts/');

%raw data location
data_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/raw_data/';

%processed data location (for saving)
proc_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/';

%Select subject to process if not specified
if nargin<1
    choice_made=0;
    sub_folders=dir([data_folder filesep 'S*']);
    fprintf('Select subject to analyse: \n')
    for s=1:size(sub_folders,1)
        fprintf(['[' int2str(s) '] ' sub_folders(s).name '\n'])
    end
    while ~choice_made
        choice=input('Enter number: ');
        if choice>size(sub_folders,1) || choice<1
            disp('Invalid choice, please try again');
        else
            sub=sub_folders(choice).name;
            disp(['Subject ' sub ' selected'])
            choice_made=1;
        end
    end
end
    
%identify session containing raw MEG data
sub_folder=[data_folder sub filesep];
subfolders=dir(sub_folder);
for i=1:length(subfolders)    
    if length(subfolders(i).name)>2
        ses_dir_files=dir([sub_folder subfolders(i).name]);
        for n=1:length(ses_dir_files)
            if ~isempty(strfind('MEG',ses_dir_files(n).name))
                raw_data=[sub_folder subfolders(i).name filesep ses_dir_files(n).name filesep];
            end
            if ~isempty(strfind('beh',ses_dir_files(n).name))
                logs=[sub_folder subfolders(i).name filesep ses_dir_files(n).name filesep];
            end
        end
    end
end

datasets=dir([raw_data '*.fif']);

%Long session can be split into multiple files; these have to be
%added in the right order.
split_cnt=1;
split_parts=zeros(1,length(datasets));
split_files=[];
single_files=[];
for d=1:length(datasets)
    has_split=strcmp([datasets(d).name(1:end-4) '-1.fif'],{datasets.name});
    if any(has_split)
        split_files{split_cnt}={datasets(d),datasets(find(has_split))};
        split_cnt=split_cnt+1;
        split_parts(d)=1;
        split_parts(find(has_split))=1;        
    end
end
if any(~split_parts)
    single_inds=find(~split_parts);
    for i=1:length(single_inds)
        single_files{i,1}=datasets(single_inds(i));
    end
end
   
%put all files in chronological order
all_dats=[split_files ; single_files];

for n=1:length(all_dats)
    if iscell(all_dats{n})
        tmp_hdr=ft_read_header([raw_data all_dats{n}{1}.name]);
    else
        tmp_hdr=ft_read_header([raw_data all_dats{n}.name]);
    end
    timestamp(n)=tmp_hdr.orig.raw.info.meas_date(1);
end
[~,order]=sort(timestamp);

datasets=[];
for n=1:length(all_dats)
    datasets{n}=all_dats{order(n)};
end

%% Import data for all datasets
for d=1:length(datasets)        
    
    %import data
    if ~iscell(datasets{d}) %single file
        disp(['Opening dataset: ' datasets{d}.name]);
        cfg=[];
        cfg.dataset = [raw_data datasets{d}.name];
        data = ft_preprocessing(cfg);
        event = ft_read_event(cfg.dataset);
    else %split file
        %load in the split dataset
        disp(['Opening split dataset: ' datasets{d}{1}.name]);
        cfg=[];
        cfg.dataset = [raw_data datasets{d}{1}.name];
        d1 = ft_preprocessing(cfg);
        cfg.dataset = [raw_data datasets{d}{2}.name];
        d2 = ft_preprocessing(cfg);
        
        %concatenate data
        disp('Concatenating data...')
        data=d1;
        data.hdr.nSamples=d1.hdr.nSamples+d2.hdr.nSamples;
        
        %concatenating time
        time_step=median(gradient(d1.time{1}));
        if (data.fsample-(1/time_step))<1e-7
            disp(['Sampling rate of ' int2str(1/time_step) 'Hz seems consistent.'])
        else
            error('Check sampling rate')
        end
        
        data.time{1}=[d1.time{1} d2.time{1}+d1.time{1}(end)+time_step];
        s1=size(d1.trial{1});
        s2=size(d2.trial{1});
        data.trial{1}=zeros(s1(1),s1(2)+s2(2));
        data.trial{1}=horzcat(d1.trial{1}, d2.trial{1});
        data.sampleinfo=[1 d1.sampleinfo(2)+d2.sampleinfo(2)];
        
        d1_file = [raw_data datasets{d}{1}.name];
        d2_file = [raw_data datasets{d}{2}.name];
        
        disp('Getting trigger events.')
        event = ft_read_event({d1_file,d2_file});
    end
        
        
    %from the data, define MEG channels
    MEGchannels=strmatch('meg',data.hdr.chantype);
    MISCchannels=strmatch('misc',data.hdr.chantype);
    EOGchannels=strmatch('BIO',data.hdr.label);
    
    %separate MAG and planar
    MEGMAGchannels=strmatch('megmag',data.hdr.chantype);
    MEGPLANARchannels=strmatch('megplanar',data.hdr.chantype);
    
    %Select MEG planar channels
    disp('Preprocessing planar sensors')
    cfg=[];
    cfg.channel = MEGPLANARchannels;
    cfg.continuous = 'yes';
    cfg.hpfilter      = 'yes';
    cfg.hpfreq        = 1;
    cfg.demean        = 'yes';
    data_MEG_planar = ft_preprocessing(cfg,data);
    
    %Select MEG magnetometer channels,demean and filter
    disp('Preprocessing magnetometers')
    cfg=[];
    cfg.channel = MEGMAGchannels;
    cfg.continuous = 'yes';    
    cfg.hpfilter      = 'yes';
    cfg.hpfreq        = 1;
    cfg.demean        = 'yes';
    data_MEG_mag = ft_preprocessing(cfg,data);
    
    %get scaling for other channels to bring into similar range
    scaling=std(data_MEG_planar.trial{1}(:));
    
    %extract misc channels (diode, eyelink)
    disp('Preprocessing miscellaneous channels')
    cfg=[];
    cfg.channel = MISCchannels;
    cfg.continuous = 'yes';
    cfg.demean        = 'yes';
    data_MISC = ft_preprocessing(cfg,data);
    
    data_MISC.label{1}='ELX';
    data_MISC.label{2}='ELY';
    data_MISC.label{3}='ELPUPIL';
    data_MISC.label{4}='DIODE';
      
    %extract EOG and filter
    disp('Preprocessing EOG')
    cfg=[];
    cfg.channel = EOGchannels;
    cfg.continuous = 'yes';
    cfg.demean        = 'yes';
    cfg.lpfilter      = 'yes'; %'no' or 'yes'  lowpass filter (default = 'no')
    cfg.hpfilter      = 'yes'; %'no' or 'yes'  highpass filter (default = 'no')
    cfg.lpfreq        = 30; %lowpass  frequency in Hz
    cfg.hpfreq        = 1; %highpass frequency in Hz
    data_EOG = ft_preprocessing(cfg,data);
    
    data_EOG.trial{1}(1,:)=data_EOG.trial{1}(1,:)*(scaling/std(data_EOG.trial{1}(1,:)));
    data_EOG.label{1}='vEOG';       
    
    %Bring it all back together
    disp('Reassembling data')
    hdr_orig=data.hdr;
    data = ft_appenddata([],data_MEG_mag,data_MEG_planar,data_MISC,data_EOG);
    
    %adjust header to appended data specifications
    hdr=hdr_orig;
    
    %correct chanunit and chantype
    hdr.chantype=repmat({'unknown'},size(data.label,1),1);
    hdr.chanunit=repmat({'unknown'},size(data.label,1),1);
    
    for c=1:length(data.label)
        ind_orig=find(strcmp(data.label{c},hdr_orig.label));
        if ~isempty(ind_orig)
            hdr.chantype{c}=hdr_orig.chantype{ind_orig};
            hdr.chanunit{c}=hdr_orig.chanunit{ind_orig};
        end
    end
    
    %adjust labels, chans and samples
    hdr.label=data.label;
    hdr.nChans=size(data.label,1);
    hdr.nSamples=size(data.trial{1},2);         
    
    %Segment data into trials for all conditions
    disp('Cutting epochs')
    ev_values=[2 4 8 16]; %Conditions 1-4
    
    %cue-locked
    cfg=[];
    cfg.event=event;
    cfg.hdr=hdr;
    cfg.trialfun                = 'ft_trialfun_general';
    cfg.trialdef.eventtype      ='Trigger';
    cfg.trialdef.eventvalue     = ev_values;
    cfg.trialdef.prestim        = 2.3; % in seconds
    cfg.trialdef.poststim       = 1; % in seconds
    cfg = ft_definetrial(cfg);
    trl=cfg.trl;
    
    %Discrimination target-locked
    cfg=[];
    cfg.event=event;
    cfg.hdr=hdr;
    cfg.trialfun                = 'ft_trialfun_general';
    cfg.trialdef.eventtype      ='Trigger';
    cfg.trialdef.eventvalue     = 32;
    cfg.trialdef.prestim        = 3; % in seconds
    cfg.trialdef.poststim       = 1; % in seconds
    cfg = ft_definetrial(cfg);
    trl_dt=cfg.trl;
    
    cfg=[];
    cfg.trl=trl;
    
    MEG_run{d}.data=ft_redefinetrial(cfg,data);
    MEG_run{d}.grad=data.grad;
    
    cfg.trl=trl_dt; 
    MEG_run{d}.data_dt=ft_redefinetrial(cfg,data);
    
    %Copy condition info from cue-locked analysis to discrimination target
    %analysis - but check if they are equal
    if size(MEG_run{d}.data.trialinfo,1)==size(MEG_run{d}.data_dt.trialinfo,1)
        MEG_run{d}.data_dt.trialinfo=MEG_run{d}.data.trialinfo;
    else
        error('Cue-locked analysis has different amount of trials that discrimination target locked. Something is wrong');
    end
    
    disp('Done for this run')

    clear data    
end

%% Sorting
%Now that we have all data, lets put everything together

disp('Getting logfiles')
logfiles=dir([logs '*.mat']);

for d=1:size(logfiles,1)
    log=load([logfiles(d).folder filesep logfiles(d).name]);
    run_nr_log=log.cfg.block;
    ConMat{run_nr_log}=log.cfg.ConMat;    
    [~,conf]=sort(log.cfg.FreqMat);
    config(run_nr_log)=conf(1);
    behavior{run_nr_log}=log.cfg.behavior;
    exp_cfg{run_nr_log}=log.cfg;
end

%Check whether the first dataset is also the first run
data=[];
for d=1:length(datasets)    
    %match logfile to dataset
    if length(datasets{d})>1
        run_nr(d)=str2num(datasets{d}{1}.name(findstr('run',lower(datasets{d}{1}.name))+3:end-4));
    else
        run_nr(d)=str2num(datasets{d}.name(findstr('run',lower(datasets{d}.name))+3:end-4));
    end
end

%combine MEG data such that run 1 is first
cfg=[];
cfg.keepsampleinfo='no'; %concatenating two files, discarding sampleinfo
data.meg=ft_appenddata(cfg,MEG_run{run_nr(1)}.data,MEG_run{run_nr(2)}.data);
data.meg_dt=ft_appenddata(cfg,MEG_run{run_nr(1)}.data_dt,MEG_run{run_nr(2)}.data_dt);

%store configuration for run 1 and 2
data.config=[ones(length(MEG_run{run_nr(1)}.data.trial),1)*config(1) ; ones(length(MEG_run{run_nr(2)}.data.trial),1)*config(2)];

%store behavior in same format
data.behavior=[behavior{1} ; behavior{2}];

%and condition matrix
data.ConMat=[ConMat{1} ; ConMat{2}];

%save experimental cfg
data.exp_cfg{1}=exp_cfg{1};
data.exp_cfg{2}=exp_cfg{2};

%save grads
data.grad{1}=MEG_run{run_nr(1)}.grad;
data.grad{2}=MEG_run{run_nr(2)}.grad;

%load eyelink files
disp('locating eyelink files..');
try
    load([proc_folder sub filesep sub '_eyelink.mat']);
    data.el=el;
    data.el_dt=el_dt;
    disp('Found ''em')
catch
    warning('No eyelink files detected! Omitting..')
end

%% Saving
fprintf(['Saving data for ' sub '...'])
if exist([proc_folder sub],'dir')>0
    save([proc_folder sub filesep sub '_proc.mat'],'data','-v7.3');
else %create dir
    disp('Creating subject directory..')
    mkdir([proc_folder sub]);
    save([proc_folder sub filesep sub '_proc.mat'],'data','-v7.3');
end
fprintf('..Done\n')

end
