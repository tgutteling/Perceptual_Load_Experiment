function [] = C1b_Group_RFT_dt()
%
% Group RFT results for discrimination target-locked analysis
% ROI created in the cue-locked analysis will be used

%use only correct trials?
correct_only = 1;

%set path
run('/rds/projects/2017/jenseno-01/Tjerk/set_path');

%folders
proc_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/';

%Check for High frequency TFRs in subject folders
sub_folders=dir([proc_folder filesep 'S*']);
cnt=1;
if correct_only
    for s=1:size(sub_folders,1)
        if exist([proc_folder sub_folders(s).name filesep sub_folders(s).name '_TFR_HF_dt_correct_only.mat'])>0
            datasets{cnt}=[proc_folder sub_folders(s).name filesep sub_folders(s).name '_TFR_HF_dt_correct_only.mat'];
            cnt=cnt+1;
        end
    end
else
    for s=1:size(sub_folders,1)
        if exist([proc_folder sub_folders(s).name filesep sub_folders(s).name '_TFR_HF_dt.mat'])>0
            datasets{cnt}=[proc_folder sub_folders(s).name filesep sub_folders(s).name '_TFR_HF_dt.mat'];
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

%channel selection
for d=1:length(data)
    MEG_sens=strmatch('MEG',data{d}.TFR.left.HF.ind{1,1}.label);    
    planars{d}=find(cellfun(@(x) ~isempty(strfind(x,'+')),data{d}.TFR.left.HF.ind{1,1}.label(MEG_sens),'UniformOutput',1));
end

%Store freq and time
freq=data{1}.TFR.left.HF.ind{1,1}.freq;
time=data{1}.TFR.left.HF.ind{1,1}.time;

%get RFT indices
ind_63=dsearchn(freq',63);
ind_70=dsearchn(freq',70);

%% Load ROIs
if correct_only
    load([proc_folder filesep 'group' filesep 'ROI_RFT_correct_only.mat'],'ROI_left','ROI_right')
else
    load([proc_folder filesep 'group' filesep 'ROI_RFT.mat'],'ROI_left','ROI_right')
end

%% Extract timecourses

for d=1:length(data)
    cur=data{d}.TFR;
    for l=1:4 %4 load conditions
        cfg=[];
                              
        %targets left
        cfg.channel=ROI_left;
        tmp=ft_selectdata(cfg,cur.left.HF.ind{1,l});
        RFT_dt{d}.ind.target_left{l,1}=squeeze(nanmean(tmp.powspctrm(:,ind_63,:),1));        
        
        tmp=ft_selectdata(cfg,cur.left.HF.ind{2,l});
        RFT_dt{d}.ind.target_left{l,2}=squeeze(nanmean(tmp.powspctrm(:,ind_70,:),1));        
        
        %targets right
        cfg.channel=ROI_right;
        tmp=ft_selectdata(cfg,cur.right.HF.ind{2,l});
        RFT_dt{d}.ind.target_right{l,1}=squeeze(nanmean(tmp.powspctrm(:,ind_63,:),1));
               
        tmp=ft_selectdata(cfg,cur.right.HF.ind{1,l});
        RFT_dt{d}.ind.target_right{l,2}=squeeze(nanmean(tmp.powspctrm(:,ind_70,:),1));
                
        %distractors left
        cfg.channel=ROI_left;
        tmp=ft_selectdata(cfg,cur.right.HF.ind{1,l});
        RFT_dt{d}.ind.dist_left{l,1}=squeeze(nanmean(tmp.powspctrm(:,ind_63,:),1));
                
        tmp=ft_selectdata(cfg,cur.right.HF.ind{2,l});
        RFT_dt{d}.ind.dist_left{l,2}=squeeze(nanmean(tmp.powspctrm(:,ind_70,:),1));
                
        %distractors right
        cfg.channel=ROI_right;
        tmp=ft_selectdata(cfg,cur.left.HF.ind{2,l});
        RFT_dt{d}.ind.dist_right{l,1}=squeeze(nanmean(tmp.powspctrm(:,ind_63,:),1));        
        
        tmp=ft_selectdata(cfg,cur.left.HF.ind{1,l});
        RFT_dt{d}.ind.dist_right{l,2}=squeeze(nanmean(tmp.powspctrm(:,ind_70,:),1));        
    end    
end

%% Saving
fprintf(['Saving data RFT group data for ' int2str(length(data)) ' subjects..'])
if ~(exist([proc_folder 'group'],'dir')>0)
    disp('Creating subject directory..')
    mkdir([proc_folder 'group']);
end

if correct_only
    save([proc_folder 'group' filesep 'ROI_RFT_dt_correct_only.mat'],'RFT_dt','RFT_dt_ami','time','ROI_left','ROI_right','ROI_left_ami','ROI_right_ami','-v7.3');
else
    save([proc_folder 'group' filesep 'ROI_RFT_dt.mat'],'RFT_dt','time','ROI_left','ROI_right','-v7.3');
end

fprintf('..Done\n')

end


