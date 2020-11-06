function [el] = A0_el_preproc(sub)
%
% Parses eyelink data from ASC file
% EDF files first have to be converted to ASC before processing
% 
% Eyelink data trials are segmeneted and sorted according to t=0 (cue
% offset). Order of trials is acconrding to  block nr (not
% counterbalancing)
%
% Current limitation in using 'prestim' time as data is segmented from
% trial onset onwards

%raw data folder
data_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/raw_data/';

%processed data location (for saving)
proc_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/';

%Subject selection if none specified
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
    
sub_folder=[data_folder sub filesep];
subfolders=dir(sub_folder);
for i=1:length(subfolders)    
    if length(subfolders(i).name)>2
        ses_dir_files=dir([sub_folder subfolders(i).name]);
        for n=1:length(ses_dir_files)
            if ~isempty(strfind('eyelink',ses_dir_files(n).name))
                el_dir=[sub_folder subfolders(i).name filesep ses_dir_files(n).name filesep];
            end
        end
    end
end
            
fls=dir([el_dir '*.asc']);       
if isempty(fls)
    error('No eyelink files found..');
end

for f=1:size(fls,1)
    el=[];
    cur=fls(f).name;
    disp(['Parsing file ' int2str(f)]);
    
    %Open eyelink ASC file
    fid= fopen([el_dir filesep cur]);
    
    %ignore preamle
    line = fgetl(fid);
    while isempty(strfind(line, 'Config'))
        line = fgetl(fid);
    end
    
    %Extract first trial info
    vals= sscanf(line, 'MSG	%i %i %i Config: %i Attention side: %i FaceNr: %i Movement direction: %i Congruency: %i');
    
    %values:
    % 1 - eyelink time stamp 
    % 2 - block
    % 3 - trial
    % 4 - Load condition
    % 5 - attention side
    % 6 - Face
    % 7 - eye movement direction
    % 8 - Congruency
    
    disp(['Parsing file ' int2str(f) '/' int2str(size(fls,1)) ' Block: ' int2str(vals(2)) ' Trial: ' int2str(vals(3))])
    
    el{1,1}.config=vals(4);
    el{1,1}.att_side=vals(5);
        
    %get to first trigger
    while isempty(strfind(line, 'Trigger'))
        line = fgetl(fid);
    end    
    
    vals= sscanf(line, 'MSG	%i %i %i %f %f Trigger: %i');
    %values:
    % 1 - eyelink time stamp
    % 2 - block
    % 3 - trial
    % 4 - block time
    % 5 - trial time
    % 6 - trigger
        
    el_t0_block=vals(1); %eyelink time zero
    %el_t0_trial=vals(1);
    block=vals(2);
    stim_PC_t0_block=vals(4); %Display PC time zero
    %stim_PC_t0_trial=vals(5);
    el{1,1}.Trigger(1,:)=[0 vals(4) vals(6)];
    
    trialend=0;
    cnt=1;
    tr_cnt=2;
    tmp=[];    
    trial=1;
    line = fgetl(fid);
    while trialend==0
        while isempty(strfind(line, 'Config')) && isempty(strfind(line, 'END')) && trialend==0
            if isempty(strfind(line, 'MSG'))
                vals=[];
                vals = sscanf(line, '%i %f %f %f %f');
                if ~isempty(vals)
                    l=length(vals);
                    if l<4
                        vals(l+1:4)=NaN; %data missing, fill in NaNs
                    end                        
                    vals(1)=vals(1)-el_t0_block; %use experiment startttime
                    try
                        tmp(cnt,:)=vals(1:4);
                    catch
                        disp('here')
                    end
                    cnt=cnt+1;
                end
                line = fgetl(fid);
            else
                if ~isempty(strfind(line, 'MSG')) && ~isempty(strfind(line, 'Trigger'))
                    vals= sscanf(line, 'MSG	%i %i %i %f %f Trigger: %i');                    
                    %el{trial}.Trigger(tr_cnt,:)=[vals(1)-el_t0_block vals(6)];
                    el{trial}.Trigger(tr_cnt,:)=[vals(1)-el_t0_block vals(4)-stim_PC_t0_block vals(6)];
                    tr_cnt=tr_cnt+1;
                end
                line = fgetl(fid);
            end
        end
        el{trial}.time=tmp(:,1);
        el{trial}.x=tmp(:,2);
        el{trial}.y=tmp(:,3);
        el{trial}.pupil=tmp(:,4);
        tmp=[];
        cnt=1;
        tr_cnt=1;        
        
        if ~isempty(strfind(line, 'END OF SESSION'))        
            trialend=1;
            break;
        end
        
        %Scan next trial info
        vals= sscanf(line, 'MSG	%i %i %i Config: %i Attention side: %i FaceNr: %i Movement direction: %i Congruency: %i');        
        trial=vals(3);
        disp(['Parsing file ' int2str(f) '/' int2str(size(fls,1)) ' Block: ' int2str(block) ' Trial: ' int2str(trial)])
        el{trial}.config=vals(4);
        el{trial}.att_side=vals(5);
        line = fgetl(fid);        
    end
    ses{block}=el;
end

%save all (raw data)
disp('Eyelink data parsed, saving...');
save([el_dir 'el_data_raw.mat'],'ses')

% meg_time=[];
% el_time=[];
% for i=1:length(el)
%     meg_time=[meg_time ; el{i}.Trigger(1,1)];
%     el_time=[el_time ; el{i}.Trigger(1,2)];
% end

%% %%%%%%%%%%%%%%%%%%%%%%%
% Create trial structure %
%%%%%%%%%%%%%%%%%%%%%%% %%

%TriggerStart = 1;   % Trial start
%TriggerFig1 = 2;    % config 1
%TriggerFig2 = 4;    % config 2
%TriggerFig3 = 8;    % config 3
%TriggerFig4 = 16;   % config 4
%TriggerMov = 32;    % Onset of stimulus eye movement
%TriggerResp = 64;   % subject response

disp('Sorting trials for cue-locked analysis')
%pre and poststimulus time, assuming 1000Hz sampling rate 
%taking cue offset (triggers 2-16) as t=0;
prestim=1300; %ms baseline, before averaging trigger
poststim=1000; %ms post stimulus onset
av_trig=2;
trig_column=size(ses{1}{1}.Trigger,2);

global_trialnr=1;
el=[];
for s=1:length(ses)
    cur=ses{s};
    %trials=size(cur,2);
    trials=sum(~cellfun('isempty',cur),2);
    for t=1:trials          
        trl=cur{t};
        if sum(ismember(trl.Trigger(av_trig,trig_column),[2,4,8,16]))>0 %incomplete are ignored
            el.config(global_trialnr)=trl.config;
            el.att_side(global_trialnr)=trl.att_side;
            el.BlockNr(global_trialnr)=s;
            el.TrialNr(global_trialnr)=t;
            if ~ismember(trl.Trigger(av_trig,trig_column),[2,4,8,16])
                warning('Inconsistent trigger');
            end
            cueindex=find(trl.time==trl.Trigger(av_trig,1));
            if size(trl.x,1)<length(cueindex-prestim:cueindex+poststim)
                error(['Eyelink recorded length too short! Data may be corrupted in session ' int2str(s) ' trial ' int2str(t)])
            else
                el.x(global_trialnr,:)=trl.x(cueindex-prestim:cueindex+poststim);
                el.y(global_trialnr,:)=trl.y(cueindex-prestim:cueindex+poststim);
                el.pupil(global_trialnr,:)=trl.pupil(cueindex-prestim:cueindex+poststim);
                global_trialnr=global_trialnr+1;
            end
        end
    end
end

%time axis is the same for all
el.time=(-prestim:poststim)/1000;

disp('Sorting trials for discrimination target-locked analysis')
%pre and poststimulus time, assuming 1000Hz sampling rate 
%taking discrimination target (trigger  32) as t=0;
prestim=2000; %ms baseline, before averaging trigger
poststim=800; %ms post stimulus onset
av_trig=3;
trig_column=size(ses{1}{1}.Trigger,2);

global_trialnr=1;
el_dt=[];
for s=1:length(ses)
    cur=ses{s};
    %trials=size(cur,2);
    trials=sum(~cellfun('isempty',cur),2);
    for t=1:trials
        trl=cur{t};
        if sum(ismember(trl.Trigger(av_trig,trig_column),32))>0 %incomplete are ignored
            el_dt.config(global_trialnr)=trl.config;
            el_dt.att_side(global_trialnr)=trl.att_side;
            el_dt.BlockNr(global_trialnr)=s;
            el_dt.TrialNr(global_trialnr)=t;
            if ~ismember(trl.Trigger(av_trig,trig_column),32)
                warning('Inconsistent trigger');
            end
            cueindex=find(trl.time==trl.Trigger(av_trig,1));
            %if size(trl.x,1)<length(cueindex-prestim:cueindex+poststim)
            if size(trl.x,1)<(cueindex+poststim)
                warning(['Eyelink recorded length too short! Data may be corrupted in session ' int2str(s) ' trial ' int2str(t)])
                trl.x(end:cueindex+poststim)=NaN; %NaN-pad
                trl.y(end:cueindex+poststim)=NaN; %NaN-pad
                trl.pupil(end:cueindex+poststim)=NaN; %NaN-pad
                el_dt.x(global_trialnr,:)=trl.x(cueindex-prestim:cueindex+poststim);
                el_dt.y(global_trialnr,:)=trl.y(cueindex-prestim:cueindex+poststim);
                el_dt.pupil(global_trialnr,:)=trl.pupil(cueindex-prestim:cueindex+poststim);
                global_trialnr=global_trialnr+1;
            else
                el_dt.x(global_trialnr,:)=trl.x(cueindex-prestim:cueindex+poststim);
                el_dt.y(global_trialnr,:)=trl.y(cueindex-prestim:cueindex+poststim);
                el_dt.pupil(global_trialnr,:)=trl.pupil(cueindex-prestim:cueindex+poststim);
                global_trialnr=global_trialnr+1;
            end
        end
    end
end

%time axis is the same for all
el_dt.time=(-prestim:poststim)/1000;

%save all (epoched data)
disp('Saving sorted data...')
if exist([proc_folder sub],'dir')>0
    save([proc_folder sub filesep  sub '_eyelink.mat'],'el','el_dt');
else %create dir
    disp('Creating subject directory..')
    mkdir([proc_folder sub]);
    save([proc_folder sub filesep  sub '_eyelink.mat'],'el','el_dt');
end

disp('Done');
        

end