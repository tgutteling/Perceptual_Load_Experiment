function [] = B1_Freq_Domain_TFR_HF(sub,av_type)
%
% Freq domain analyses - TFRs - High frequencies
%

correct_only=1; %use only correct trials

%set path
run('/rds/projects/2017/jenseno-01/Tjerk/set_path');

%folders
proc_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/';

%subject/condition selection when not specified
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

%Load Processed data
switch av_type
    case 1
        fprintf(['Loading ' proc_folder sub filesep sub '_all_clean.mat...'])
        load([proc_folder sub filesep sub '_all_clean.mat']);disp('Done')
        
    case 2
        fprintf(['Loading ' proc_folder sub filesep sub '_all_clean_dt.mat...'])
        load([proc_folder sub filesep sub '_all_clean_dt.mat']);disp('Done')
end


%% HF TFR

disp('[1] Computing high frequency TFR')

switch av_type
    case 1 %Cue-locked
        %[1] high frequencies
        cfg              = [];
        cfg.foi          = 50:0.5:100;
        cfg.pad          = 4;
        cfg.toi          = -2.3:0.01:1;
        cfg.output       = 'pow';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.t_ftimwin    = repmat(0.5,1,length(cfg.foi));%10./cfg.foi;%
        cfg.channel      = strmatch('MEG',data{1}.left{1}.meg.label);
        
        for c=1:2 %both configurations
            for l=1:4 %4 load conditions
                
                disp(['Calculating TFR for configuration ' int2str(c) ', attention left, load condition ' int2str(l)])
                if correct_only
                    cfg.trials = find(data{c}.left{l}.behavior(:,10)>0);
                    beh.left{c,l}=data{c}.left{l}.behavior(cfg.trials,:);
                else
                    beh.left{c,l}=data{c}.left{l}.behavior;
                end
                TFR.left.HF.ind{c,l} = ft_freqanalysis(cfg,data{c}.left{l}.meg);
                TFR.left.HF.ind{c,l}.grad = data{c}.grad;
                
                if correct_only
                    cfg.trials = find(data{c}.right{l}.behavior(:,10)>0);
                    beh.right{c,l}=data{c}.right{l}.behavior(cfg.trials,:);
                else
                    beh.right{c,l}=data{c}.right{l}.behavior;
                end
                
                disp(['Calculating TFR for configuration ' int2str(c) ', attention right, load condition ' int2str(l)])
                TFR.right.HF.ind{c,l} = ft_freqanalysis(cfg,data{c}.right{l}.meg);
                TFR.right.HF.ind{c,l}.grad = data{c}.grad;                                
            end
        end
        
    case 2 %Discrimination target-locked
        %[1] high frequencies
        cfg              = [];
        cfg.foi          = 50:0.5:100;
        cfg.pad          = 4;%'nextpow2';
        cfg.toi          = -2.3:0.01:1;
        cfg.output       = 'pow';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.t_ftimwin    = repmat(0.5,1,length(cfg.foi));
        cfg.channel      = strmatch('MEG',data{1}.left{1}.meg_dt.label);
        
        for c=1:2 %both configurations
            for l=1:4 %4 load conditions
                
                disp(['Calculating TFR for configuration ' int2str(c) ', attention left, load condition ' int2str(l)])
                if correct_only
                    cfg.trials = find(data{c}.left{l}.behavior(:,10)>0);
                    beh.left{c,l}=data{c}.left{l}.behavior(cfg.trials,:);
                else
                    beh.left{c,l}=data{c}.left{l}.behavior;
                end
                TFR.left.HF.ind{c,l} = ft_freqanalysis(cfg,data{c}.left{l}.meg_dt);
                TFR.left.HF.ind{c,l}.grad = data{c}.grad;
                
                if correct_only
                    cfg.trials = find(data{c}.right{l}.behavior(:,10)>0);
                    beh.right{c,l}=data{c}.right{l}.behavior(cfg.trials,:);
                else
                    beh.right{c,l}=data{c}.right{l}.behavior;
                end
                
                disp(['Calculating TFR for configuration ' int2str(c) ', attention right, load condition ' int2str(l)])
                TFR.right.HF.ind{c,l} = ft_freqanalysis(cfg,data{c}.right{l}.meg_dt);
                TFR.right.HF.ind{c,l}.grad = data{c}.grad;
            end
        end
end
        

%combine planar data
disp('Creating combined planar data')
for c=1:2 %both configurations
    for l=1:4 %4 load conditions
        TFR.left.HF.ind{c,l}=ft_combineplanar([],TFR.left.HF.ind{c,l});
        TFR.right.HF.ind{c,l}=ft_combineplanar([],TFR.right.HF.ind{c,l});
    end
end


%% RFT timecourse

disp('[1] Computing TFR power timecourse')

%[1] high frequencies
cfg              = [];
cfg.foi          = 63;
cfg.pad          = 4;
cfg.toi          = -2.3:0.01:1;
cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.keeptrials   = 'yes';
cfg.t_ftimwin    = repmat(0.5,1,length(cfg.foi));
if av_type==1
    cfg.channel      = strmatch('MEG',data{1}.left{1}.meg.label);
else
    cfg.channel      = strmatch('MEG',data{1}.left{1}.meg_dt.label);
end

for c=1:2 %both configurations
    for l=1:4 %4 load conditions
        
        %left
        cfg.foi=63;
        disp(['Calculating 63Hz RFT trials for configuration ' int2str(c) ', attention left, load condition ' int2str(l)])
        if correct_only
            cfg.trials = find(data{c}.left{l}.behavior(:,10)>0);
        end
        if av_type==1
            TFR_trials.left.f1{c,l} = ft_freqanalysis(cfg,data{c}.left{l}.meg);
        else
            TFR_trials.left.f1{c,l} = ft_freqanalysis(cfg,data{c}.left{l}.meg_dt);
        end
        TFR_trials.left.f1{c,l}.grad = data{c}.grad;
        
        disp(['Calculating 70Hz RFT trials for configuration ' int2str(c) ', attention left, load condition ' int2str(l)])
        cfg.foi = 70;
        if av_type==1
            TFR_trials.left.f2{c,l} = ft_freqanalysis(cfg,data{c}.left{l}.meg);
        else
            TFR_trials.left.f2{c,l} = ft_freqanalysis(cfg,data{c}.left{l}.meg_dt);
        end
        TFR_trials.left.f2{c,l}.grad = data{c}.grad;
        
        %right
        if correct_only
            cfg.trials = find(data{c}.right{l}.behavior(:,10)>0);
        end
        
        cfg.foi=63;
        disp(['Calculating 63Hz RFT trials for configuration ' int2str(c) ', attention right, load condition ' int2str(l)])
        if av_type==1
            TFR_trials.right.f1{c,l} = ft_freqanalysis(cfg,data{c}.right{l}.meg);
        else
            TFR_trials.right.f1{c,l} = ft_freqanalysis(cfg,data{c}.right{l}.meg_dt);
        end
        TFR_trials.right.f1{c,l}.grad = data{c}.grad;
        
        cfg.foi = 70;
        disp(['Calculating 70Hz RFT trials for configuration ' int2str(c) ', attention right, load condition ' int2str(l)])
        if av_type==1
            TFR_trials.right.f2{c,l} = ft_freqanalysis(cfg,data{c}.right{l}.meg);
        else
            TFR_trials.right.f2{c,l} = ft_freqanalysis(cfg,data{c}.right{l}.meg_dt);
        end
        TFR_trials.right.f2{c,l}.grad = data{c}.grad;
    end
end

%combine planar data
disp('Creating combined planar data')
for c=1:2 %both configurations
    for l=1:4 %4 load conditions
        TFR_trials.left.f1{c,l}=ft_combineplanar([],TFR_trials.left.f1{c,l});
        TFR_trials.right.f1{c,l}=ft_combineplanar([],TFR_trials.right.f1{c,l});
        
        TFR_trials.left.f2{c,l}=ft_combineplanar([],TFR_trials.left.f2{c,l});
        TFR_trials.right.f2{c,l}=ft_combineplanar([],TFR_trials.right.f2{c,l});
    end
end

%% Save data
fprintf(['Saving data for ' sub '...'])
switch av_type
    case 1
        if correct_only
            save([proc_folder sub filesep sub '_TFR_HF_correct_only.mat'],'TFR','TFR_trials','beh','-v7');
        else
            save([proc_folder sub filesep sub '_TFR_HF.mat'],'TFR','TFR_trials','beh','-v7');
        end
    case 2
        if correct_only
            save([proc_folder sub filesep sub '_TFR_HF_dt_correct_only.mat'],'TFR','TFR_trials','beh','-v7');
        else
            save([proc_folder sub filesep sub '_TFR_HF_dt.mat'],'TFR','TFR_trials','beh','-v7');
        end
end
fprintf('..Done\n')

end


