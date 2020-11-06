function [] = C2a_Group_Alpha()
%
% Group level Alpha TFR and timecourse extraction for the cue-locked data. ROI is calculated from the
% low frequency data
%

%use only correct trials?
correct_only = 1;

%the rejected subject should not be part of the ROI forming
bad_subs=[17,18,20,24,29];

%set path
run('/rds/projects/2017/jenseno-01/Tjerk/set_path');

%folders
proc_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/';

%Check for Low frequency TFRs in subject folders
sub_folders=dir([proc_folder filesep 'S*']);
cnt=1;
if correct_only
    for s=1:size(sub_folders,1)
        if exist([proc_folder sub_folders(s).name filesep sub_folders(s).name '_TFR_LF_correct_only.mat'])>0
            datasets{cnt}=[proc_folder sub_folders(s).name filesep sub_folders(s).name '_TFR_LF_correct_only.mat'];
            cnt=cnt+1;
        end
    end
else
    for s=1:size(sub_folders,1)
        if exist([proc_folder sub_folders(s).name filesep sub_folders(s).name '_TFR_LF.mat'])>0
            datasets{cnt}=[proc_folder sub_folders(s).name filesep sub_folders(s).name '_TFR_LF.mat'];
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
for d=1:length(datasets)
    MEG_sens=strmatch('MEG',data{d}.TFR.left.LF.ind{1,1}.label);    
    planars{d}=find(cellfun(@(x) ~isempty(strfind(x,'+')),data{d}.TFR.left.LF.ind{1,1}.label(MEG_sens),'UniformOutput',1));
end

%% Create TFR GA subtraction for plotting
for d=1:length(datasets)
    cfg=[];
    cfg.channel=planars{d};
    att_left{d}=ft_freqgrandaverage(cfg,data{d}.TFR.left.LF.ind{:,:});
    att_right{d}=ft_freqgrandaverage(cfg,data{d}.TFR.right.LF.ind{:,:});
end

%create GA
att_left_GA=ft_freqgrandaverage([],att_left{:});
att_right_GA=ft_freqgrandaverage([],att_right{:});

%remove cfg for storage
att_left_GA.cfg=[];
att_right_GA.cfg=[];

%save
disp('Saving Grand average TFR')
save([proc_folder 'group' filesep 'LF_TFR_GA_correct_only_figplot.mat'],'att_left_GA','att_right_GA')


%% ROI selection
load neuromag306cmb.mat

%Store time
time=data{1}.TFR.left.LF.ind{1,1}.time;

%Nr of sensors to select
n_best=4;

%Create average per subject
for d=1:size(data,2)
    
    %collapse all to attention left and right
    left = ft_appendfreq([],data{d}.TFR.left.LF.ind{:,:});
    right = ft_appendfreq([],data{d}.TFR.right.LF.ind{:,:});
    
    %get list of sensors available for all subjects
    common_labels=intersect(data{1}.TFR.left.LF.ind{1,1}.label(planars{1}), data{2}.TFR.left.LF.ind{1,1}.label(planars{2}));
    for i=3:size(data,2)
        common_labels=intersect(common_labels, data{i}.TFR.left.LF.ind{1,1}.label(planars{i}));
    end        
    
    %select planars and alpha range
    cfg=[];
    cfg.channel=common_labels;
    cfg.frequency=[8 13];
    cfg.latency=[0 1];
    cfg.avgoverrpt  = 'yes';
    cfg.avgoverfreq = 'yes';
    lp{d}=ft_selectdata(cfg,left);
    rp{d}=ft_selectdata(cfg,right);   
end

%% Create group average and select ROI
good_subs=setdiff(1:35,bad_subs);
lp_av=ft_freqgrandaverage([],lp{good_subs});
rp_av=ft_freqgrandaverage([],rp{good_subs});

cfg=[];
cfg.keepindividual='yes';
lp_all_av=ft_freqgrandaverage(cfg,lp{good_subs});
rp_all_av=ft_freqgrandaverage(cfg,rp{good_subs});
lp_all_av.powspctrm=squeeze(nanmean(lp_all_av.powspctrm,4));
rp_all_av.powspctrm=squeeze(nanmean(rp_all_av.powspctrm,4));
p_labels=lp_av.label;

rp_chanpower=squeeze(nanmean(rp_av.powspctrm,3));
lp_chanpower=squeeze(nanmean(lp_av.powspctrm,3));

%now get sensors with best modulation
[~,ind]=sort(rp_chanpower-lp_chanpower);
best_planar_left_names = flipud(p_labels(ind(end-(n_best-1):end)));
[~,ind]=sort(lp_chanpower-rp_chanpower);
best_planar_right_names  = flipud(p_labels(ind(end-(n_best-1):end)));

%% Create plot data to check sensor selection
sub=[];
sub.label=lp_all_av.label;
sub.avg=(lp_all_av.powspctrm-rp_all_av.powspctrm)';
sub.time=1;
sub.dimord='subj_chan_time';

plot_data.att_left=sub;

sub.avg=(rp_all_av.powspctrm-lp_all_av.powspctrm)';
plot_data.att_right=sub;
plot_data.roi_left=best_planar_left_names;
plot_data.roi_right=best_planar_right_names;

%% Extract timecourses from ROIs
for d=1:length(data)
    cur=data{d}.TFR;
    for l=1:4 %4 load conditions
        cfg=[];
        disp(['Subject ' int2str(d) ', Load condition ' int2str(l)])
        
        ROI_left=best_planar_left_names;
        ROI_right=best_planar_right_names;
        
        %targets left
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
        
        %distractors right
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
    save([proc_folder 'group' filesep 'ROI_alpha_correct_only.mat'],'alpha','time','ROI_left','ROI_right','plot_data','-v7.3');
else
    save([proc_folder 'group' filesep 'ROI_alpha_' int2str(type) '.mat'],'alpha','time','ROI_left','ROI_right','plot_data','-v7.3');
end

fprintf('..Done\n')

end


