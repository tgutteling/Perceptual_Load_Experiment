function [] = B2_Freq_Domain_TFR_LF(sub,av_type)
%
% Freq domain analyses - TFRs - Low frequencies
% -Usage 
%       B4_Freq_Domain_TFR_LF(sub,av_type)
%
%    where
%        sub      - Subject to analyse, string (e.g. 'S01')
%        av_type  - Averaging type, cue-locked (1) or discrimination
%                   target-locked (2)

%use only correct trials?
correct_only=1;

%set path
run('/rds/projects/2017/jenseno-01/Tjerk/set_path');

%folders
proc_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/';

%slect type if not specified
if nargin<2
    choice_made=0;
    fprintf('Select averaging type: \n')
    fprintf('[1] Cue-locked\n')
    fprintf('[2] Discrimination target-locked\n')
    while ~choice_made
        choice=input('Enter number: ');
        if choice>2 || choice<1
            disp('Invalid choice, please try again');
        else
            av_type=choice;
            switch choice
                case 1
                    disp('Cue-locked selected')
                case 2
                    disp('Discrimination target-locked selected')
            end
            choice_made=1;
        end
    end
end

%select subject if not specified
if nargin<1
    choice_made=0;
    sub_folders=dir([proc_folder filesep 'S*']);
    fprintf('Select subject to analyse: \n')
    switch av_type
        case 1
            cnt=1;
            for s=1:size(sub_folders,1)
                if exist([proc_folder sub_folders(s).name filesep sub_folders(s).name '_all_clean.mat'])>0
                    fprintf(['[' int2str(cnt) '] ' sub_folders(s).name '\n'])
                    sel(cnt)=s;
                    cnt=cnt+1;
                end
            end
        case 2
            cnt=1;
            for s=1:size(sub_folders,1)
                if exist([proc_folder sub_folders(s).name filesep sub_folders(s).name '_all_clean_dt.mat'])>0
                    fprintf(['[' int2str(cnt) '] ' sub_folders(s).name '\n'])
                    sel(cnt)=s;
                    cnt=cnt+1;
                end
            end
    end
    while ~choice_made
        choice=input('Enter number: ');
        if choice>(cnt-1) || choice<1
            disp('Invalid choice, please try again');
        else
            sub=sub_folders(sel(choice)).name;
            disp(['Subject ' sub ' selected'])
            choice_made=1;
        end
    end
end

switch av_type
    case 1
        fprintf(['Loading ' proc_folder sub filesep sub '_all_clean.mat...'])
        load([proc_folder sub filesep sub '_all_clean.mat']);disp('Done')
    case 2
        fprintf(['Loading ' proc_folder sub filesep sub '_all_clean_dt.mat...'])
        load([proc_folder sub filesep sub '_all_clean_dt.mat']);disp('Done')
end


%% LF TFR

disp('[1] Computing low frequency TFR')

%[1] low frequencies
cfg              = [];
cfg.foi          = 2:30;
cfg.pad          = 'nextpow2';
cfg.toi          = -2.3:0.01:1;
cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.t_ftimwin    = 3./cfg.foi;
if av_type==1
    cfg.channel      = strmatch('MEG',data{1}.left{1}.meg.label);
else
    cfg.channel      = strmatch('MEG',data{1}.left{1}.meg_dt.label);
end

beh=[];
for c=1:2 %both configurations
    for l=1:4 %4 load conditions
        disp(['Calculating TFR for configuration ' int2str(c) ', attention left, load condition ' int2str(l)])
        
        %attention left
        if correct_only
            cfg.trials = find(data{c}.left{l}.behavior(:,10)>0);
            beh.left{c,l}=data{c}.left{l}.behavior(cfg.trials,:);
        else
            beh.left{c,l}=data{c}.left{l}.behavior;
        end
        if av_type==1
            TFR.left.LF.ind{c,l} = ft_freqanalysis(cfg,data{c}.left{l}.meg);
        else
            TFR.left.LF.ind{c,l} = ft_freqanalysis(cfg,data{c}.left{l}.meg_dt);
        end
        TFR.left.LF.ind{c,l}.grad = data{c}.grad;
        
        %attention right
        if correct_only
            cfg.trials = find(data{c}.right{l}.behavior(:,10)>0);
            beh.right{c,l}=data{c}.right{l}.behavior(cfg.trials,:);
        else
            beh.right{c,l}=data{c}.right{l}.behavior;
        end    
        
        disp(['Calculating TFR for configuration ' int2str(c) ', attention right, load condition ' int2str(l)])
        if av_type==1
            TFR.right.LF.ind{c,l} = ft_freqanalysis(cfg,data{c}.right{l}.meg);
        else
            TFR.right.LF.ind{c,l} = ft_freqanalysis(cfg,data{c}.right{l}.meg_dt);
        end
        TFR.right.LF.ind{c,l}.grad = data{c}.grad;
    end
end

%combine planar data
disp('Creating combined planar data')
for c=1:2 %both configurations
    for l=1:4 %4 load conditions
        TFR.left.LF.ind{c,l}=ft_combineplanar([],TFR.left.LF.ind{c,l});        
        TFR.right.LF.ind{c,l}=ft_combineplanar([],TFR.right.LF.ind{c,l});        
    end
end

%% Just for alpha timecourses
%[1] low frequencies
cfg              = [];
cfg.foi          = 8:13;
cfg.pad          = 'nextpow2';
cfg.toi          = -2.3:0.01:1;
cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.t_ftimwin    = 3./cfg.foi;
cfg.keeptrials   = 'yes';
if av_type==1
    cfg.channel      = strmatch('MEG',data{1}.left{1}.meg.label);
else
    cfg.channel      = strmatch('MEG',data{1}.left{1}.meg_dt.label);
end

for c=1:2 %both configurations
    for l=1:4 %4 load conditions
        disp(['Calculating TFR trials for configuration ' int2str(c) ', attention left, load condition ' int2str(l)])
        
        %attention left
        if correct_only
            cfg.trials = find(data{c}.left{l}.behavior(:,10)>0);            
        end
        if av_type==1
            TFR_trials.left.LF{c,l} = ft_freqanalysis(cfg,data{c}.left{l}.meg);
        else
            TFR_trials.left.LF{c,l} = ft_freqanalysis(cfg,data{c}.left{l}.meg_dt);
        end
        TFR_trials.left.LF{c,l}.grad = data{c}.grad;
        
        %attention right
        if correct_only
            cfg.trials = find(data{c}.right{l}.behavior(:,10)>0);
        end    
        
        disp(['Calculating TFR trials for configuration ' int2str(c) ', attention right, load condition ' int2str(l)])
        if av_type==1
            TFR_trials.right.LF{c,l} = ft_freqanalysis(cfg,data{c}.right{l}.meg);
        else
            TFR_trials.right.LF{c,l} = ft_freqanalysis(cfg,data{c}.right{l}.meg_dt);
        end
        TFR_trials.right.LF{c,l}.grad = data{c}.grad;
    end
end

%average over frequency
cfg=[];
cfg.avgoverfreq = 'yes';
for c=1:2 %both configurations
    for l=1:4 %4 load conditions
        TFR_trials.left.LF{c,l}=ft_selectdata(cfg,TFR_trials.left.LF{c,l});        
        TFR_trials.right.LF{c,l}=ft_selectdata(cfg,TFR_trials.right.LF{c,l});        
    end
end
        
%combine planar data
disp('Creating combined planar data')
for c=1:2 %both configurations
    for l=1:4 %4 load conditions
        TFR_trials.left.LF{c,l}=ft_combineplanar([],TFR_trials.left.LF{c,l});        
        TFR_trials.right.LF{c,l}=ft_combineplanar([],TFR_trials.right.LF{c,l});        
    end
end

    
%% Save data    
fprintf(['Saving data for ' sub '...'])
switch av_type
    case 1
        if correct_only
            save([proc_folder sub filesep sub '_TFR_LF_correct_only.mat'],'TFR','TFR_trials','beh','-v7.3');
        else
            save([proc_folder sub filesep sub '_TFR_LF.mat'],'TFR','TFR_trials','beh','-v7.3');
        end
    case 2
        if correct_only
            save([proc_folder sub filesep sub '_TFR_LF_dt_correct_only.mat'],'TFR','TFR_trials','beh','-v7.3');
        else
            save([proc_folder sub filesep sub '_TFR_LF_dt.mat'],'TFR','TFR_trials','beh','-v7.3');
        end
end
fprintf('..Done\n')

end


