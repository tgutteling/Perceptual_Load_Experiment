function [] = C1a_Group_RFT()
%
% Group level RFT TFR and timecourse extraction. ROI is calculated from the
% high frequency data
%

%use only correct trials?
correct_only = 1;

%set path
run('/rds/projects/2017/jenseno-01/Tjerk/set_path');

%folders
proc_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/';

%the rejected subjects should not be part of the ROI forming
bad_subs=[17,18,20,24,29];

%Check for high frequency TFRs in subject folders
sub_folders=dir([proc_folder filesep 'S*']);
cnt=1;
if correct_only
    for s=1:size(sub_folders,1)
        if exist([proc_folder sub_folders(s).name filesep sub_folders(s).name '_TFR_HF_correct_only.mat'])>0
            datasets{cnt}=[proc_folder sub_folders(s).name filesep sub_folders(s).name '_TFR_HF_correct_only.mat'];
            sub_nr(cnt)=str2num(sub_folders(s).name(2:end));
            cnt=cnt+1;
        end
    end
else
    for s=1:size(sub_folders,1)
        if exist([proc_folder sub_folders(s).name filesep sub_folders(s).name '_TFR_HF.mat'])>0
            datasets{cnt}=[proc_folder sub_folders(s).name filesep sub_folders(s).name '_TFR_HF.mat'];
            sub_nr(cnt)=str2num(sub_folders(s).name(2:end));
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
    MEG_sens=strmatch('MEG',data{d}.TFR.left.HF.ind{1,1}.label);    
    planars{d}=find(cellfun(@(x) ~isempty(strfind(x,'+')),data{d}.TFR.left.HF.ind{1,1}.label(MEG_sens),'UniformOutput',1));
end

%% Create TFR GA subtraction for plotting
for d=1:length(datasets)
    cfg=[];
    cfg.channel=planars{d};
    att63_left{d}=ft_freqgrandaverage(cfg,data{d}.TFR.left.HF.ind{1,:});
    att70_left{d}=ft_freqgrandaverage(cfg,data{d}.TFR.left.HF.ind{2,:});
    att63_right{d}=ft_freqgrandaverage(cfg,data{d}.TFR.right.HF.ind{2,:});
    att70_right{d}=ft_freqgrandaverage(cfg,data{d}.TFR.right.HF.ind{1,:});
end

%create GA
att63_left_GA=ft_freqgrandaverage([],att63_left{:});
att63_right_GA=ft_freqgrandaverage([],att63_right{:});
att70_left_GA=ft_freqgrandaverage([],att70_left{:});
att70_right_GA=ft_freqgrandaverage([],att70_right{:});

%remove cfg as it is irrelevant and takes up a lot of space
att63_left_GA.cfg=[];
att63_right_GA.cfg=[];
att70_left_GA.cfg=[];
att70_right_GA.cfg=[];

%save
disp('Saving Grand average TFR')
save([proc_folder 'group' filesep 'HF_TFR_GA_correct_only_figplot.mat'],'att63_left_GA','att63_right_GA','att70_left_GA','att70_right_GA')


%% ROI creation
% Here we create the roi based on the maximal group modulation

%Store freq and time
freq=data{1}.TFR.left.HF.ind{1,1}.freq;
time=data{1}.TFR.left.HF.ind{1,1}.time;

%get RFT indices
ind_63=dsearchn(freq',63);
ind_70=dsearchn(freq',70);

%get inds for baseline
bl=[-1 -0.4];
ind0=find(time>=bl(1),1,'first');
ind1=find(time>=bl(2),1,'first');

%find last non-nan datapoint
t_end=data{1}.TFR.left.HF.ev{1,1}.time(find(~isnan(squeeze(data{1}.TFR.left.HF.ev{1,1}.powspctrm(1,1,:))),1,'last'));

%get list of sensors available for all subjects
common_labels=intersect(data{1}.TFR.left.HF.ind{1,1}.label(planars{1}), data{2}.TFR.left.HF.ind{1,1}.label(planars{2}));
for d=3:size(data,2)
    common_labels=intersect(common_labels, data{d}.TFR.left.HF.ind{1,1}.label(planars{d}));
end

%extract data per direction and target/distractor
%select only planar grads
for d=1:size(data,2)
    %collapse over load
    c1_left = ft_appendfreq([],data{d}.TFR.left.HF.ev{1,:});
    c1_right = ft_appendfreq([],data{d}.TFR.right.HF.ev{1,:});
    c2_left = ft_appendfreq([],data{d}.TFR.left.HF.ev{2,:});
    c2_right = ft_appendfreq([],data{d}.TFR.right.HF.ev{2,:});
    
    cfg=[];
    cfg.latency=[0 t_end];
    cfg.frequency=[63 63];
    cfg.avgovertime='yes';
    cfg.avgoverfreq = 'yes';
    cfg.avgoverrpt='yes';
    cfg.channel=common_labels;%planars{d};
    all{d}.left_target{1}=ft_selectdata(cfg,c1_left);
    all{d}.right_target{1}=ft_selectdata(cfg,c2_right);
    all{d}.left_dist{1}=ft_selectdata(cfg,c1_right);
    all{d}.right_dist{1}=ft_selectdata(cfg,c2_left);
    
    cfg.frequency=[70 70];
    all{d}.left_target{2}=ft_selectdata(cfg,c2_left);
    all{d}.right_target{2}=ft_selectdata(cfg,c1_right);
    all{d}.left_dist{2}=ft_selectdata(cfg,c2_right);
    all{d}.right_dist{2}=ft_selectdata(cfg,c1_left);
end

%collapso to single matrix
for d=1:length(all)
    left_target_63(d,:)=all{d}.left_target{1}.powspctrm;
    right_target_63(d,:)=all{d}.right_target{1}.powspctrm;
    left_dist_63(d,:)=all{d}.left_dist{1}.powspctrm;
    right_dist_63(d,:)=all{d}.right_dist{1}.powspctrm;
    
    left_target_70(d,:)=all{d}.left_target{2}.powspctrm;
    right_target_70(d,:)=all{d}.right_target{2}.powspctrm;
    left_dist_70(d,:)=all{d}.left_dist{2}.powspctrm;
    right_dist_70(d,:)=all{d}.right_dist{2}.powspctrm;
end


%Select best sensors
n_best=4;

%remove bad subs
left_target_63(bad_subs,:)=[];
right_target_63(bad_subs,:)=[];
left_dist_63(bad_subs,:)=[];
right_dist_63(bad_subs,:)=[];
left_target_70(bad_subs,:)=[];
right_target_70(bad_subs,:)=[];
left_dist_70(bad_subs,:)=[];
right_dist_70(bad_subs,:)=[];

%Demean and calculate group mean
left_target_63=left_target_63./mean(left_target_63,2);m_lt63=mean(left_target_63,1);
right_target_63=right_target_63./mean(right_target_63,2);m_rt63=mean(right_target_63,1);
left_dist_63=left_dist_63./mean(left_dist_63,2);m_ld63=mean(left_dist_63,1);
right_dist_63=right_dist_63./mean(right_dist_63,2);m_rd63=mean(right_dist_63,1);

left_target_70=left_target_70./mean(left_target_70,2);m_lt70=mean(left_target_70,1);
right_target_70=right_target_70./mean(right_target_70,2);m_rt70=mean(right_target_70,1);
left_dist_70=left_dist_70./mean(left_dist_70,2);m_ld70=mean(left_dist_70,1);
right_dist_70=right_dist_70./mean(right_dist_70,2);m_rd70=mean(right_dist_70,1);

m_left_target=mean([m_lt63 ; m_lt70]);
m_right_target=mean([m_rt63 ; m_rt70]);
m_left_dist=mean([m_ld63 ; m_ld70]);
m_right_dist=mean([m_rd63 ; m_rd70]);

all_left_target=(left_target_63 + left_target_70)/2;
all_right_target=(right_target_63 + right_target_70)/2;
all_left_dist=(left_dist_63 + left_dist_70)/2;
all_right_dist=(right_dist_63 + right_dist_70)/2;

%now get sensors with best modulation
[~,ind]=sort(m_left_target-m_left_dist);
best_planar_left_names = flipud(common_labels(ind(end-(n_best-1):end)));
[~,ind]=sort(m_right_target-m_right_dist);
best_planar_right_names  = flipud(common_labels(ind(end-(n_best-1):end)));

%Create sensor level average for plotting
sub=[];
sub.label=common_labels;
sub.avg=all_left_target';
sub.time=1;
sub.dimord='subj_chan_time';

plot_data.left_target=sub;
sub.avg=all_right_target';
plot_data.right_target=sub;
sub.avg=all_left_dist';
plot_data.left_dist=sub;
sub.avg=all_right_dist';
plot_data.right_dist=sub;
sub.avg=ami_left';
plot_data.ami_left=sub;
sub.avg=ami_right';
plot_data.ami_right=sub;

plot_data.labels=common_labels;
plot_data.roi_left=best_planar_left_names;
plot_data.roi_right=best_planar_right_names;
plot_data.roi_left_ami=best_planar_left_names_ami;
plot_data.roi_right_ami=best_planar_right_names_ami;


%% Extract timecourses
for d=1:length(data)
    cur=data{d}.TFR;
    for l=1:4 %4 load conditions
        cfg=[];
        
        ROI_left=best_planar_left_names;
        ROI_right=best_planar_right_names;                
        
        %targets
        cfg.channel=ROI_left;
        tmp=ft_selectdata(cfg,cur.left.HF.ind{1,l});
        RFT{d}.ind.target_left{l,1}=squeeze(nanmean(tmp.powspctrm(:,ind_63,:),1));
        
        tmp=ft_selectdata(cfg,cur.left.HF.ind{2,l});
        RFT{d}.ind.target_left{l,2}=squeeze(nanmean(tmp.powspctrm(:,ind_70,:),1));
        
        %targets right
        cfg.channel=ROI_right;
        tmp=ft_selectdata(cfg,cur.right.HF.ind{2,l});
        RFT{d}.ind.target_right{l,1}=squeeze(nanmean(tmp.powspctrm(:,ind_63,:),1));
        
        tmp=ft_selectdata(cfg,cur.right.HF.ind{1,l});
        RFT{d}.ind.target_right{l,2}=squeeze(nanmean(tmp.powspctrm(:,ind_70,:),1));
        
        %distractors left
        cfg.channel=ROI_left;
        tmp=ft_selectdata(cfg,cur.right.HF.ind{1,l});
        RFT{d}.ind.dist_left{l,1}=squeeze(nanmean(tmp.powspctrm(:,ind_63,:),1));
        
        tmp=ft_selectdata(cfg,cur.right.HF.ind{2,l});
        RFT{d}.ind.dist_left{l,2}=squeeze(nanmean(tmp.powspctrm(:,ind_70,:),1));
        
        %distractors right
        cfg.channel=ROI_right;
        tmp=ft_selectdata(cfg,cur.left.HF.ind{2,l});
        RFT{d}.ind.dist_right{l,1}=squeeze(nanmean(tmp.powspctrm(:,ind_63,:),1));
        
        tmp=ft_selectdata(cfg,cur.left.HF.ind{1,l});
        RFT{d}.ind.dist_right{l,2}=squeeze(nanmean(tmp.powspctrm(:,ind_70,:),1));                
    end
end

%% Saving
fprintf(['Saving data RFT group data for ' int2str(length(data)) ' subjects..'])
if ~(exist([proc_folder 'group'],'dir')>0)
    disp('Creating subject directory..')
    mkdir([proc_folder 'group']);
end

if correct_only
    save([proc_folder 'group' filesep 'ROI_RFT_correct_only.mat'],'RFT','time','ROI_left','ROI_right','plot_data','-v7.3');
else
    save([proc_folder 'group' filesep 'ROI_RFT.mat'],'RFT','time','ROI_left','ROI_right','plot_data','-v7.3');
end
fprintf('..Done\n')

end


