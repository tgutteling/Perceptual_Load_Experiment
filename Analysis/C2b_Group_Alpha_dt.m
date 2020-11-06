function [] = C2_Group_Alpha_dt()
%
% Group level Alpha TFR and timecourse extraction for the discrimination target-locked data. ROI is calculated from the
% low frequency data

%use only correct trials?
correct_only = 1;

%set path
run('/rds/projects/2017/jenseno-01/Tjerk/set_path');

%folders
proc_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/';

%Check for Low frequency TFRs in subject folders
sub_folders=dir([proc_folder filesep 'S*']);
cnt=1;
if correct_only
    for s=1:size(sub_folders,1)
        if exist([proc_folder sub_folders(s).name filesep sub_folders(s).name '_TFR_LF_dt_correct_only.mat'])>0
            datasets{cnt}=[proc_folder sub_folders(s).name filesep sub_folders(s).name '_TFR_LF_dt_correct_only.mat'];
            cnt=cnt+1;
        end
    end
else
    for s=1:size(sub_folders,1)
        if exist([proc_folder sub_folders(s).name filesep sub_folders(s).name '_TFR_LF_dt.mat'])>0
            datasets{cnt}=[proc_folder sub_folders(s).name filesep sub_folders(s).name '_TFR_LF_dt.mat'];
            cnt=cnt+1;
        end
    end
end

%load data
disp(['Found ' int2str(length(datasets)) ' processed datasets'])
for d=1:length(datasets)
    disp(['Loading ' datasets{d}])
    data{d}=load(datasets{d});
end

%% ROI selection, defined by cue-locked data
if correct_only
    load([proc_folder filesep 'group' filesep 'ROI_alpha_correct_only.mat'],'ROI_left','ROI_right');
else
    load([proc_folder filesep 'group' filesep 'ROI_alpha.mat'],'ROI_left','ROI_right');
end

%Store time
time=data{1}.TFR.left.LF.ind{1,1}.time;

%% Get timecourses from ROI

for d=1:length(data)
    cur=data{d}.TFR;
    for l=1:4 %4 load conditions
        cfg=[];
        disp(['Subject ' int2str(d) ', Load condition ' int2str(l)])
               
        %target left
        tmp=ft_appendfreq([],cur.left.LF.ind{:,l});
        cfg.channel=ROI_left;
        cfg.frequency =[8 13];
        cfg.avgoverchan = 'yes';
        cfg.avgoverfreq='yes';
        tmp=ft_selectdata(cfg,tmp);
        alpha{d}.target_left{l}=squeeze(mean(tmp.powspctrm,1));                
        
        %targets right
        tmp=ft_appendfreq([],cur.right.LF.ind{:,l});
        cfg.channel=ROI_right;
        tmp=ft_selectdata(cfg,tmp);
        alpha{d}.target_right{l}=squeeze(mean(tmp.powspctrm,1));
                       
        %distractors left
        tmp=ft_appendfreq([],cur.right.LF.ind{:,l});
        cfg.channel=ROI_left;
        cfg.frequency =[8 13];
        cfg.avgoverchan = 'yes';
        cfg.avgoverfreq='yes';
        tmp=ft_selectdata(cfg,tmp);
        alpha{d}.dist_left{l}=squeeze(mean(tmp.powspctrm,1));
               
        %targets right
        tmp=ft_appendfreq([],cur.left.LF.ind{:,l});
        cfg.channel=ROI_right;
        tmp=ft_selectdata(cfg,tmp);
        alpha{d}.dist_right{l}=squeeze(mean(tmp.powspctrm,1));                
    end
end

%% Saving
fprintf(['Saving data alpha group data for ' int2str(length(data)) ' subjects..'])
if ~exist([proc_folder 'group'],'dir')>0 
    disp('Creating subject directory..')
    mkdir([proc_folder 'group']);    
end
if correct_only
    save([proc_folder 'group' filesep 'ROI_alpha_dt_correct_only.mat'],'alpha','time','ROI_left','ROI_right','-v7.3');
else
    save([proc_folder 'group' filesep 'ROI_alpha_dt.mat'],'alpha','time','ROI_left','ROI_right','-v7.3');
end

fprintf('..Done\n')

end


