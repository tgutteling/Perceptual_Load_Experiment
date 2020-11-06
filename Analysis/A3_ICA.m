function [] = A3_ICA(sub,av_type)
%
% Removal of bad sensors and trials is applied and ICA components are
% calculated and stored for further inspection
% av_type: 1 - cue-locked 2 - Discrimination target-locked

%% Set up folders

%set path
run('/rds/projects/2017/jenseno-01/Tjerk/set_path');

%proc data location
data_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/';

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

if nargin<2
    av_type=1; %Cue-locked
end

%load data
disp(['Loading ' data_folder sub filesep sub '_proc.mat'])
load([data_folder sub filesep sub '_proc.mat']);

%get MEG channels
switch av_type
    case 1 %%Cue-locked
        disp('Selecting MEG channels..')
        MEGchannels=strmatch('MEG',data.meg.label);
        cfg=[];
        cfg.channel=MEGchannels;
        all=ft_selectdata(cfg,data.meg);
        data=rmfield(data,'meg_dt');
    case 2 %Discrimination target-locked
        disp('Selecting MEG channels..')
        MEGchannels=strmatch('MEG',data.meg_dt.label);
        cfg=[];
        cfg.channel=MEGchannels;
        all=ft_selectdata(cfg,data.meg_dt);
        data=rmfield(data,'meg');
end

%% Load results of visual inspection
try
    switch av_type
        case 1
            load([data_folder sub filesep sub '_preICA_bad.mat']);
        case 2
            load([data_folder sub filesep sub '_preICA_bad_dt.mat']);
    end
catch
    error('No pre ICA cleanup performed, aborting..')
end


%% Remove the worst from data
trials=1:length(all.trial);
trials(bad_trials)=[];
cfg=[];
cfg.trials=trials;
cfg.channel={'all'};
for b=1:length(bad_chans)
    cfg.channel{b+1}=['-' bad_chans{b}];
end
all=ft_selectdata(cfg,all);

%do the same for all the synced data
data.config(bad_trials)=[];
data.behavior(bad_trials,:)=[];
data.ConMat(bad_trials,:)=[];

%eyelink
switch av_type
    case 1
        data.el.config(bad_trials)=[];
        data.el.att_side(bad_trials)=[];
        data.el.BlockNr(bad_trials)=[];
        data.el.TrialNr(bad_trials)=[];
        data.el.x(bad_trials,:)=[];
        data.el.y(bad_trials,:)=[];
        data.el.pupil(bad_trials,:)=[];
    case 2
        data.el_dt.config(bad_trials)=[];
        data.el_dt.att_side(bad_trials)=[];
        data.el_dt.BlockNr(bad_trials)=[];
        data.el_dt.TrialNr(bad_trials)=[];
        data.el_dt.x(bad_trials,:)=[];
        data.el_dt.y(bad_trials,:)=[];
        data.el_dt.pupil(bad_trials,:)=[];
end

%save aux channels
switch av_type
    case 1        
        AUXchannels=[strmatch('EL',data.meg.label) ; strmatch('DIODE',data.meg.label) ; strmatch('vEOG',data.meg.label)];
        cfg=[];
        cfg.channel=AUXchannels;
        trials=1:length(data.meg.trial);
        trials(bad_trials)=[];
        cfg.trials=trials;
        aux=ft_selectdata(cfg,data.meg);
    case 2
        AUXchannels=[strmatch('EL',data.meg_dt.label) ; strmatch('DIODE',data.meg_dt.label) ; strmatch('vEOG',data.meg_dt.label)];
        cfg=[];
        cfg.channel=AUXchannels;
        trials=1:length(data.meg_dt.trial);
        trials(bad_trials)=[];
        cfg.trials=trials;
        aux=ft_selectdata(cfg,data.meg_dt);
end

%separate MEG from AUX channels for recompose
if av_type==1    
    data.meg=all;
else
    data.meg_dt=all;
end
data.aux=aux;

%% ICA
disp('Starting ICA..')
cfg        = [];  
cfg.method = 'runica';  
cfg.runica.maxsteps = 100;  
comp = ft_componentanalysis(cfg,all);   

data.comp=comp;

%save comps for later reference/offline analysis
disp('Done! Saving components')
switch av_type
    case 1
        save([data_folder sub filesep sub '_comp.mat'],'data','-v7.3');
    case 2
        save([data_folder sub filesep sub '_comp_dt.mat'],'data','-v7.3');
end
disp('Also Done!')

end


