function [] = D2_Paper_figures()
%
% Here we only plot all the relevant figures as reported in the paper
%

%folders
group_folder='Z:\Tjerk\Load2\proc_data\group\';
addpath Z:\Tjerk\Scripts_general\Violinplot\
addpath Z:\Tjerk\Scripts_general\
addpath Z:\Tjerk\Load2\Layouts\
addpath Z:\fieldtrip-latest\


%% subject rejection
bad_subs=[17,18,20,24,29];

%Figure 1 is the paradigm overview

%% Figure 2 - Behavioral results
% Violin plot of the Reaction times per load, and the distractor effect
f2=figure('Name','Fig 2: Behaviour');
load([group_folder 'Group_RT.mat']);

%remove bad subs
RT_all(bad_subs,:)=[];

%fancy violin plot - RT
beh_plot=[RT_all(:,1) RT_all(:,3) RT_all(:,2) RT_all(:,4)];

subplot(1,3,[1 2])
h=violinplot(beh_plot,{'noisy distractors','salient distractors','noisy distractors','salient distractors'},'ShowMean',true);
h(1,1).ViolinColor=[0.2 0.2 0.8];
h(1,1).EdgeColor=[0.2 0.2 0.8];
h(1,1).ViolinPlot.LineStyle='--';
h(1,1).ViolinPlot.LineWidth=2;
h(1,2).ViolinColor=[0.2 0.2 0.8];
h(1,2).EdgeColor=[0.2 0.2 0.8];
h(1,2).ViolinPlot.LineWidth=2;
h(1,3).ViolinColor=[0.8 0.2 0.2];
h(1,3).EdgeColor=[0.8 0.2 0.2];
h(1,3).ViolinPlot.LineWidth=2;
h(1,3).ViolinPlot.LineStyle='--';
h(1,4).ViolinColor=[0.8 0.2 0.2];
h(1,4).EdgeColor=[0.8 0.2 0.2];
h(1,4).ViolinPlot.LineWidth=2;
ylabel('Reaction time (ms)');ylim([300 650])
xlabel('Low target load                               High target load')
text(-0.1,650,'A','Fontsize',20)

subplot(1,3,3)
beh_plot=[RT_all(:,3)-RT_all(:,1) RT_all(:,4)-RT_all(:,2)];

h=violinplot(beh_plot,{'Low load','high load'},'ShowMean',true);
h(1,1).ViolinColor=[0.2 0.2 0.8];
h(1,1).EdgeColor=[0.2 0.2 0.8];
h(1,1).ViolinPlot.LineWidth=2;
h(1,2).ViolinColor=[0.8 0.2 0.2];
h(1,2).EdgeColor=[0.8 0.2 0.2];
h(1,2).ViolinPlot.LineWidth=2;
ylabel('Distractor interference (ms)')
text(0.92,55,'*','Fontsize',30)
text(-0.05,60,'B','Fontsize',20);hline(0,'k');

%% Figure 3 - sensor selection & attentional effect (Alpha + RFT)
f3=figure('Name','Fig 3: Sensor selection');

%Subplot positions
sub_pos(1,:)=[0.42 0.8 0.2 0.13];%alpha topo
sub_pos(2,:)=[0.07 0.8 0.28 0.13];%left TFR alpha
sub_pos(3,:)=[0.69 0.8 0.28 0.13];%right TFR alpha
sub_pos(4,:)=[0.07 0.57 0.43 0.18]; %Alpha cue-locked timecourse
sub_pos(5,:)=[0.54 0.57 0.43 0.18]; %Alpha DT-locked timecourse

sub_pos(6,:)=[0.42 0.39 0.2 0.13]; %RFT topo att left
sub_pos(7,:)=[0.42 0.24 0.2 0.13];%RFT topo att right
sub_pos(8,:)=[0.07 0.33 0.28 0.13]; %RFT TFR left
sub_pos(9,:)=[0.69 0.33 0.28 0.13]; %RFT TFR right
sub_pos(10,:)=[0.07 0.055 0.43 0.18]; %RFT Cue-locked timecourse
sub_pos(11,:)=[0.54 0.055 0.43 0.18]; %RFT DT-locked timecourse

%Load plot data from ROI creation
tmp=load([group_folder 'ROI_alpha_4_correct_only.mat']);
load('neuromag306cmb.mat');

topo_alpha=tmp.plot_data.att_left;
topo_alpha.avg=tmp.plot_data.att_left.avg-tmp.plot_data.att_right.avg;
topo_alpha.avg=mean(topo_alpha.avg,2);
ROI_left_indx=[];
for r=1:length(tmp.plot_data.roi_left)
    ROI_left_indx(r)=strmatch(tmp.plot_data.roi_left(r),tmp.plot_data.att_left.label);
end

ROI_right_indx=[];
for r=1:length(tmp.plot_data.roi_right)
    ROI_right_indx(r)=strmatch(tmp.plot_data.roi_right(r),tmp.plot_data.att_right.label);
end

alpha.ROI_left=tmp.plot_data.roi_left;
alpha.ROI_right=tmp.plot_data.roi_right;

%Plot Alpha ROI topo
subplot('position',sub_pos(1,:))
cfg=[];
cfg.layout = lay;
cfg.gridscale = 200;
cfg.style='straight';
cfg.zlim=[-2.5e-24 2.5e-24];
cfg.comment='no';
cfg.markercolor =[0.3 0.3 0.3];

ft_topoplotER(cfg,topo_alpha);hold on;
ft_plot_lay(lay,'chanindx',[ROI_left_indx],'label','no','box','no','pointsymbol','.','pointcolor',[0 0 0],'pointsize',4)
ft_plot_lay(lay,'chanindx',[ROI_left_indx],'label','no','box','no','pointsymbol','o','pointcolor',[0 0 0],'pointsize',6)
ft_plot_lay(lay,'chanindx',[ROI_right_indx],'label','no','box','no','pointsymbol','.','pointcolor',[0 0 0],'pointsize',4)
ft_plot_lay(lay,'chanindx',[ROI_right_indx],'label','no','box','no','pointsymbol','o','pointcolor',[0 0 0],'pointsize',6)

title('Alpha (left-right)');%colorbar('NorthOutside');

%Plot Alpha TFR
load([group_folder 'LF_TFR_GA_correct_only_figplot.mat']);

LF_sub=att_left_GA;
LF_sub.powspctrm=att_left_GA.powspctrm-att_right_GA.powspctrm;
LF_sub.time=LF_sub.time+0.35;

cfg=[];
cfg.xlim=[-0.25 1];
cfg.ylim=[4 30];
cfg.zlim=[-2e-24 2e-24];
cfg.channel=alpha.ROI_right;

%TFR alpha left
subplot('position',sub_pos(2,:))
ft_singleplotTFR(cfg,LF_sub);
title('Left Hemisphere');vline(0,'k');text(0.02,5.5,'CUE onset');colorbar off
xlabel('Time (s)');ylabel('Frequency (Hz)');
text(-0.55,30,'A','Fontsize',20)

%TFR alpha right
subplot('position',sub_pos(3,:))
cfg.channel=alpha.ROI_left;
ft_singleplotTFR(cfg,LF_sub);
title('Right Hemisphere');vline(0,'k');text(0.02,5.5,'CUE onset');colorbar off
xlabel('Time (s)');ylabel('Frequency (Hz)');

%Alpha Cue-Locked Timecourse
subplot('position',sub_pos(4,:))
alpha_target_all=[];
alpha_dist_all=[];
load([group_folder 'Timecourses.mat']);
ami_alpha_all=(timecourses.alpha.ami{1}+timecourses.alpha.ami{2}+timecourses.alpha.ami{3}+timecourses.alpha.ami{4})/4;
s=size(timecourses.alpha.target);
for i=1:s(1)
    alpha_target_all(i,:)=mean([timecourses.alpha.target{i,:}],2)';
    alpha_dist_all(i,:)=mean([timecourses.alpha.dist{i,:}],2)';
end

%remove NaNs
t1=find(~isnan(ami_alpha_all(1,:)),1,'first');
t2=find(~isnan(ami_alpha_all(1,:)),1,'last');
ami_alpha_all=ami_alpha_all(:,t1:t2);
alpha_target_all=alpha_target_all(:,t1:t2);
alpha_dist_all=alpha_dist_all(:,t1:t2);
m_target_alpha=mean(alpha_target_all,1)*100;
m_dist_alpha=mean(alpha_dist_all,1)*100;
se_target_alpha=std(alpha_target_all)/sqrt(size(alpha_target_all,1))*100;
se_dist_alpha=std(alpha_dist_all)/sqrt(size(alpha_dist_all,1))*100;

cl_time=timecourses.time(t1:t2)+0.35;xlim([-0.25 cl_time(end)]);
hline(100,'color','k','linewidth',1);hold on
h1=fill([cl_time fliplr(cl_time)],[m_target_alpha-se_target_alpha fliplr(m_target_alpha+se_target_alpha)],[0.2 0.8 0.2],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([cl_time fliplr(cl_time)],[m_dist_alpha-se_dist_alpha fliplr(m_dist_alpha+se_dist_alpha)],'k','Edgealpha',0);set(h2,'FaceAlpha',0.25);hold on
h4=plot(cl_time,m_target_alpha,'color',[0.2 0.8 0.2],'linewidth',2);
h5=plot(cl_time,m_dist_alpha,'color','k','linewidth',2);
xlim([-0.25 cl_time(end)]);
ylim([90 135]);
vline(0,'color','k','linestyle','--','linewidth',2);
text(0.03,93,'CUE onset');
xlabel('Time(s)');ylabel('Alpha power (%baseline)')
legend([h4,h5],{'Target','Distractor'},'location','northeast')
text(-0.47,135,'B','Fontsize',20)

%Alpha Discrimination target-locked timecourse
subplot('position',sub_pos(5,:))
ami_dt_alpha_all=(timecourses.alpha.ami_dt{1}+timecourses.alpha.ami_dt{2}+timecourses.alpha.ami_dt{3}+timecourses.alpha.ami_dt{4})/4;
alpha_target_dt_all=[];
alpha_dist_dt_all=[];
s=size(timecourses.alpha.target);
for i=1:s(1)
    alpha_target_dt_all(i,:)=mean([timecourses.alpha.target_dt{i,:}],2)';
    alpha_dist_dt_all(i,:)=mean([timecourses.alpha.dist_dt{i,:}],2)';
end

t1=find(~isnan(ami_dt_alpha_all(1,:)),1,'first');
t2=find(~isnan(ami_dt_alpha_all(1,:)),1,'last');
alpha_target_dt_all=alpha_target_dt_all(:,t1:t2);
alpha_dist_dt_all=alpha_dist_dt_all(:,t1:t2);
m_target_dt_alpha=mean(alpha_target_dt_all,1)*100;
m_dist_dt_alpha=mean(alpha_dist_dt_all,1)*100;
se_target_dt_alpha=std(alpha_target_dt_all)/sqrt(size(alpha_target_dt_all,1))*100;
se_dist_dt_alpha=std(alpha_dist_dt_all)/sqrt(size(alpha_dist_dt_all,1))*100;
dt_time=timecourses.time(t1:t2);
xlim([-1 0.25]);
hline(100,'color','k','linewidth',1); hold on
h1=fill([dt_time fliplr(dt_time)],[m_target_dt_alpha-se_target_dt_alpha fliplr(m_target_dt_alpha+se_target_dt_alpha)],[0.2 0.8 0.2],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([dt_time fliplr(dt_time)],[m_dist_dt_alpha-se_dist_dt_alpha fliplr(m_dist_dt_alpha+se_dist_dt_alpha)],'k','Edgealpha',0);set(h2,'FaceAlpha',0.25);hold on
h4=plot(dt_time,m_target_dt_alpha,'color',[0.2 0.8 0.2],'linewidth',2);
h5=plot(dt_time,m_dist_dt_alpha,'color','k','linewidth',2);
xlim([-1 0.25]);
ylim([90 135]);
vline(0,'color','k','linestyle','--','linewidth',2);
text(-0.42,93,'Discrimination target');
xlabel('Time(s)');ax=gca;ax.YTickLabel={};

%RFT topos
tmp=load([group_folder 'ROI_RFT_type_4_correct_only.mat']);

plot_left=tmp.plot_data.left_target;
plot_left.avg=tmp.plot_data.left_target.avg-tmp.plot_data.left_dist.avg;
plot_left.avg=mean(plot_left.avg,2);
plot_right=tmp.plot_data.right_target;
plot_right.avg=tmp.plot_data.right_target.avg-tmp.plot_data.right_dist.avg;
plot_right.avg=mean(plot_right.avg,2);

ROI_left_indx=[];
for r=1:length(tmp.plot_data.roi_left)
    ROI_left_indx(r)=strmatch(tmp.plot_data.roi_left(r),tmp.plot_data.left_target.label);
end

ROI_right_indx=[];
for r=1:length(tmp.plot_data.roi_right)
    ROI_right_indx(r)=strmatch(tmp.plot_data.roi_right(r),tmp.plot_data.right_target.label);
end

RFT.ROI_left=tmp.plot_data.roi_left;
RFT.ROI_right=tmp.plot_data.roi_right;

cfg=[];
cfg.layout = lay;
cfg.gridscale = 200;
cfg.style='straight';
cfg.zlim=[-0.35 0.35];
cfg.comment='no';
cfg.markercolor =[0.3 0.3 0.3];

%RFT topo att left
subplot('position',sub_pos(6,:))
ft_topoplotER(cfg,plot_left);hold on;
ft_plot_lay(lay,'chanindx',[ROI_left_indx],'label','no','box','no','pointsymbol','.','pointcolor',[0 0 0],'pointsize',4)
ft_plot_lay(lay,'chanindx',[ROI_left_indx],'label','no','box','no','pointsymbol','o','pointcolor',[0 0 0],'pointsize',6)
title('RFT (attend left)');

%RFT topo att right
subplot('position',sub_pos(7,:))
ft_topoplotER(cfg,plot_right);
ft_plot_lay(lay,'chanindx',[ROI_right_indx],'label','no','box','no','pointsymbol','.','pointcolor',[0 0 0],'pointsize',4)
ft_plot_lay(lay,'chanindx',[ROI_right_indx],'label','no','box','no','pointsymbol','o','pointcolor',[0 0 0],'pointsize',6)
title('RFT (attend right)');

%RFT TFRs
load([group_folder 'HF_TFR_GA_correct_only_figplot.mat']);

%create subtractions
TFR_time=att63_left_GA.time+0.35;

%baseline
t1=dsearchn(TFR_time',-0.5);
t2=dsearchn(TFR_time',0);

att63_left_GA.powspctrm=(att63_left_GA.powspctrm-(mean(mean(att63_left_GA.powspctrm(:,:,t1:t2),3),2)))./(mean(mean(att63_left_GA.powspctrm(:,:,t1:t2),3),2));
att63_right_GA.powspctrm=(att63_right_GA.powspctrm-(mean(mean(att63_right_GA.powspctrm(:,:,t1:t2),3),2)))./(mean(mean(att63_right_GA.powspctrm(:,:,t1:t2),3),2));
att70_left_GA.powspctrm=(att70_left_GA.powspctrm-(mean(mean(att70_left_GA.powspctrm(:,:,t1:t2),3),2)))./(mean(mean(att70_left_GA.powspctrm(:,:,t1:t2),3),2));
att70_right_GA.powspctrm=(att70_right_GA.powspctrm-(mean(mean(att70_right_GA.powspctrm(:,:,t1:t2),3),2)))./(mean(mean(att70_right_GA.powspctrm(:,:,t1:t2),3),2));

sub_left=att63_left_GA;
sub_left.powspctrm=att63_left_GA.powspctrm-att70_left_GA.powspctrm;
sub_left.time=sub_left.time+0.35;

sub_right=att70_right_GA;
sub_right.powspctrm=att70_right_GA.powspctrm-att63_right_GA.powspctrm;
sub_right.time=sub_right.time+0.35;

cfg=[];
cfg.xlim=[-0.25 1];
cfg.ylim=[50 90];
cfg.zlim=[-0.3 0.3];
cfg.baseline=[-0.5 0];
cfg.baselinetype='absolute';
cfg.channel=[RFT.ROI_left ; RFT.ROI_right];

%RFT TFR left
subplot('position',sub_pos(8,:))
ft_singleplotTFR(cfg,sub_left);
title('Attend 63Hz');vline(0,'k');text(0.02,53,'CUE onset');colorbar off
xlabel('Time (s)');ylabel('Frequency (Hz)');
text(-0.55,90,'C','Fontsize',20)

%RFT TFR right
subplot('position',sub_pos(9,:))
cfg.channel=[RFT.ROI_left ; RFT.ROI_right];
ft_singleplotTFR(cfg,sub_right);
title('Attend 70Hz');vline(0,'k');text(0.02,53,'CUE onset');colorbar off
xlabel('Time (s)');ylabel('Frequency (Hz)');

%RFT Cue-locked timecourse
subplot('position',sub_pos(10,:))

RFT_target_all=[];
RFT_dist_all=[];
RFT_target_dt_all=[];
RFT_dist_dt_all=[];

s=size(timecourses.RFT.target);
for i=1:s(1)
    RFT_target_all(i,:)=mean([timecourses.RFT.target{i,:}],2)';
    RFT_dist_all(i,:)=mean([timecourses.RFT.dist{i,:}],2)';
    RFT_target_dt_all(i,:)=mean([timecourses.RFT.target_dt{i,:}],2)';
    RFT_dist_dt_all(i,:)=mean([timecourses.RFT.dist_dt{i,:}],2)';
end

%remove NaNs
t1=find(~isnan(mean(RFT_target_all)),1,'first');
t2=find(~isnan(mean(RFT_target_all)),1,'last');

RFT_target_all=RFT_target_all(:,t1:t2);
RFT_dist_all=RFT_dist_all(:,t1:t2);
m_target_RFT=mean(RFT_target_all,1)*100;
m_dist_RFT=mean(RFT_dist_all,1)*100;
se_target_RFT=std(RFT_target_all)/sqrt(size(RFT_target_all,1))*100;
se_dist_RFT=std(RFT_dist_all)/sqrt(size(RFT_dist_all,1))*100;

xlim([-0.25 cl_time(end)]);hline(100,'color','k','linewidth',1);hold on

%RFT Cue-locked timecourse
cl_time=timecourses.time(t1:t2)+0.35;
h1=fill([cl_time fliplr(cl_time)],[m_target_RFT-se_target_RFT fliplr(m_target_RFT+se_target_RFT)],[0.2 0.8 0.2],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([cl_time fliplr(cl_time)],[m_dist_RFT-se_dist_RFT fliplr(m_dist_RFT+se_dist_RFT)],'k','Edgealpha',0);set(h2,'FaceAlpha',0.25);hold on
h4=plot(cl_time,m_target_RFT,'color',[0.2 0.8 0.2],'linewidth',2);
h5=plot(cl_time,m_dist_RFT,'color','k','linewidth',2);
xlim([-0.25 cl_time(end)]);
ylim([45 110]);
vline(0,'color','k','linestyle','--','linewidth',1);
text(0.03,50,'CUE onset')%;hline(0,'color','k','linewidth',2);
xlabel('Time(s)');ylabel('RFT power (%baseline)')
legend([h4,h5],{'Target','Distractor'},'location','northeast')
text(-0.46,110,'D','Fontsize',20)

%RFT Discrimination target-locked timecourse
t1=find(~isnan(mean(RFT_target_dt_all)),1,'first');
t2=find(~isnan(mean(RFT_target_dt_all)),1,'last');
RFT_target_dt_all=RFT_target_dt_all(:,t1:t2);
RFT_dist_dt_all=RFT_dist_dt_all(:,t1:t2);
m_target_dt_RFT=mean(RFT_target_dt_all,1)*100;
m_dist_dt_RFT=mean(RFT_dist_dt_all,1)*100;
se_target_dt_RFT=std(RFT_target_dt_all)/sqrt(size(RFT_target_dt_all,1))*100;
se_dist_dt_RFT=std(RFT_dist_dt_all)/sqrt(size(RFT_dist_dt_all,1))*100;
dt_time=timecourses.dt_time(t1:t2);

subplot('position',sub_pos(11,:))
xlim([-1 0.2]); hline(100,'color','k','linewidth',1); hold on
h1=fill([dt_time fliplr(dt_time)],[m_target_dt_RFT-se_target_dt_RFT fliplr(m_target_dt_RFT+se_target_dt_RFT)],[0.2 0.8 0.2],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([dt_time fliplr(dt_time)],[m_dist_dt_RFT-se_dist_dt_RFT fliplr(m_dist_dt_RFT+se_dist_dt_RFT)],'k','Edgealpha',0);set(h2,'FaceAlpha',0.25);hold on
h4=plot(dt_time,m_target_dt_RFT,'color',[0.2 0.8 0.2],'linewidth',2);
h5=plot(dt_time,m_dist_dt_RFT,'color','k','linewidth',2);
xlim([-1 0.2]);
ylim([45 110]);
vline(0,'color','k','linestyle','--','linewidth',2);
text(-0.4,50,'Discrimination target');
xlabel('Time(s)');
ax=gca;ax.YTickLabel={};


%% Figure 4 - Alpha Load effects
f4=figure('Name','Fig 4: Alpha load effects');

%target alpha time window
cl_tw=[0.5 1.35];
dt_tw=[-1 0];

cl_tw_ind=[dsearchn(0.35+timecourses.time',cl_tw(1)) dsearchn(0.35+timecourses.time',cl_tw(2))];
dt_tw_ind=[dsearchn(timecourses.dt_time',dt_tw(1)) dsearchn(timecourses.dt_time',dt_tw(2))];

for i=1:4
    tmp=[timecourses.alpha.dist{:,i}]';
    dist_cl(:,i)=nanmean(tmp(:,cl_tw_ind(1):cl_tw_ind(2)),2);
    
    tmp=[timecourses.alpha.target{:,i}]';
    target_cl(:,i)=nanmean(tmp(:,cl_tw_ind(1):cl_tw_ind(2)),2);
    
    tmp=[timecourses.alpha.dist_dt{:,i}]';
    dist_dt(:,i)=nanmean(tmp(:,dt_tw_ind(1):dt_tw_ind(2)),2);
    
    tmp=[timecourses.alpha.target_dt{:,i}]';
    target_dt(:,i)=nanmean(tmp(:,dt_tw_ind(1):dt_tw_ind(2)),2);
end

%Alpha ANOVA results TARGET
subplot(3,6,1:3)
beh_plot=[target_cl(:,1) target_cl(:,3) target_cl(:,2) target_cl(:,4)];

h=violinplot(-100+beh_plot*100,{'Dist. Noisy','Dist. Salient','Dist. Noisy','Dist. Salient'},'ShowMean',true);
h(1,1).ViolinColor=[0.2 0.2 0.8];
h(1,1).EdgeColor=[0.2 0.2 0.8];
h(1,1).ViolinPlot.LineStyle='--';
h(1,1).ViolinPlot.LineWidth=2;
h(1,2).ViolinColor=[0.2 0.2 0.8];
h(1,2).EdgeColor=[0.2 0.2 0.8];
h(1,2).ViolinPlot.LineWidth=2;
h(1,3).ViolinColor=[0.8 0.2 0.2];
h(1,3).EdgeColor=[0.8 0.2 0.2];
h(1,3).ViolinPlot.LineWidth=2;
h(1,3).ViolinPlot.LineStyle='--';
h(1,4).ViolinColor=[0.8 0.2 0.2];
h(1,4).EdgeColor=[0.8 0.2 0.2];
h(1,4).ViolinPlot.LineWidth=2;
ylabel('Alpha power (rel. baseline, %change)');ylim([-40 100]);hline(0,'k');
xlabel('Low target load                 high target load')
title('Target Alpha')
text(-1.2,100,'A','Fontsize',20)

%Alpha ANOVA results DISTRACTOR
subplot(3,6,4:6)
beh_plot=[dist_cl(:,1) dist_cl(:,3) dist_cl(:,2) dist_cl(:,4)];

h=violinplot(-100+beh_plot*100,{'Dist. Noisy','Dist. Salient','Dist. Noisy','Dist. Salient'},'ShowMean',true);
h(1,1).ViolinColor=[0.2 0.2 0.8];
h(1,1).EdgeColor=[0.2 0.2 0.8];
h(1,1).ViolinPlot.LineStyle='--';
h(1,1).ViolinPlot.LineWidth=2;
h(1,2).ViolinColor=[0.2 0.2 0.8];
h(1,2).EdgeColor=[0.2 0.2 0.8];
h(1,2).ViolinPlot.LineWidth=2;
h(1,3).ViolinColor=[0.8 0.2 0.2];
h(1,3).EdgeColor=[0.8 0.2 0.2];
h(1,3).ViolinPlot.LineWidth=2;
h(1,3).ViolinPlot.LineStyle='--';
h(1,4).ViolinColor=[0.8 0.2 0.2];
h(1,4).EdgeColor=[0.8 0.2 0.2];
h(1,4).ViolinPlot.LineWidth=2;
ylim([-40 100]);hline(0,'k');ax=gca;ax.YTickLabel={};
xlabel('Low target load                 high target load')
title('Distractor Alpha')

%add ANOVA significance
line([2 4],[85 85],'color','k','linewidth',2) %low vs high load salient dist
text(3,88,'*','Fontsize',16)
line([3 4],[70 70],'color','k','linewidth',2) %high load noisy vs salient dist
text(3.5,73,'*','Fontsize',16)

%Alpha distractor time course
Alpha_dist=timecourses.alpha.dist;
Alpha_dt_dist=timecourses.alpha.dist_dt;

%Get MEANs and SEs of the timecourses
Nsubs=size(Alpha_dist,1);

%truncate to non-nan (aka numbers)
cl_t1=find(~isnan(Alpha_dist{1,1}),1,'first');
cl_t2=find(~isnan(Alpha_dist{1,1}),1,'last');

cl_time=timecourses.time(cl_t1:cl_t2)+.35; %cue onset = 0

dt_t1=find(~isnan(Alpha_dt_dist{1,1}),1,'first');
dt_t2=find(~isnan(Alpha_dt_dist{1,1}),1,'last');

dt_time=timecourses.dt_time(dt_t1:dt_t2);

selection=1:Nsubs;

%get mean and SE per load condition
for l=1:4
    
    %alpha
    m_Alpha_distractor{l}=mean([Alpha_dist{selection,l}],2);m_Alpha_distractor{l}=m_Alpha_distractor{l}(cl_t1:cl_t2)*100;
    se_Alpha_distractor{l}=std([Alpha_dist{selection,l}]')./sqrt(Nsubs);se_Alpha_distractor{l}=se_Alpha_distractor{l}(cl_t1:cl_t2)*100;
    
    %Alpha DT
    m_Alpha_dt_distractor{l}=mean([Alpha_dt_dist{selection,l}],2);m_Alpha_dt_distractor{l}=m_Alpha_dt_distractor{l}(dt_t1:dt_t2)*100;
    se_Alpha_dt_distractor{l}=std([Alpha_dt_dist{selection,l}]')./sqrt(Nsubs);se_Alpha_dt_distractor{l}=se_Alpha_dt_distractor{l}(dt_t1:dt_t2)*100;
end

%Alpha timecourses
%cue-locked
subplot(3,6,7:9)
h1=fill([cl_time fliplr(cl_time)],[m_Alpha_distractor{3}'-se_Alpha_distractor{3} fliplr(m_Alpha_distractor{3}'+se_Alpha_distractor{3})],[0 0 1],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([cl_time fliplr(cl_time)],[m_Alpha_distractor{4}'-se_Alpha_distractor{4} fliplr(m_Alpha_distractor{4}'+se_Alpha_distractor{4})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(cl_time,m_Alpha_distractor{3},'color',[0 0 1],'linewidth',2);
h4=plot(cl_time,m_Alpha_distractor{4},'color',[1 0 0],'linewidth',2);
xlim([-0.25 cl_time(end)]);ylim([90 140]);vline(0,'color','k','linestyle','--','linewidth',2);
text(0.05,94,'CUE onset');title('Distractor Alpha (salient distractors)')
xlabel('Time(s)');ylabel('Alpha power (%Baseline)')
legend([h3,h4],{'low target load','high target load'},'location','northeast')
text(-0.57,140,'B','Fontsize',20)

%discrimination target locked
subplot(3,6,10:12)
h1=fill([dt_time fliplr(dt_time)],[m_Alpha_dt_distractor{3}'-se_Alpha_dt_distractor{3} fliplr(m_Alpha_dt_distractor{3}'+se_Alpha_dt_distractor{3})],[0 0 1],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([dt_time fliplr(dt_time)],[m_Alpha_dt_distractor{4}'-se_Alpha_dt_distractor{4} fliplr(m_Alpha_dt_distractor{4}'+se_Alpha_dt_distractor{4})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(dt_time,m_Alpha_dt_distractor{3},'color',[0 0 1],'linewidth',2);
h4=plot(dt_time,m_Alpha_dt_distractor{4},'color',[1 0 0],'linewidth',2);
xlim([-cl_time(end) 0.25]);ylim([90 140]);vline(0,'color','k','linestyle','--','linewidth',2);
text(-0.55,94,'Discrimination target');
xlabel('Time(s)');ax=gca;ax.YTickLabel={};

%Add empty space for source loc
subplot(3,6,13:16)
text(-0.15,1,'C','Fontsize',20)
box off;axis off

%correlation plot
subplot(3,6,17:18)
text(-0.15,1,'D','Fontsize',20)
t1=dsearchn(cl_time',0.5)+cl_t1-1;
t2=dsearchn(cl_time',1.35)+cl_t1-1;
Alpha_LMI_salient_tw=mean(timecourses.alpha.LMI_dist_salient(:,t1:t2),2);
DE_low=RT_all(:,3)-RT_all(:,1);
scatter(Alpha_LMI_salient_tw,DE_low);
h=lsline;set(h,'color','r');hline(0,'k');vline(0,'k');
xlabel('Distractor Alpha - Load modulation index');ylabel('Beh. low load dist. interference (ms)')
[rho_sp,p_sp]=corr(Alpha_LMI_salient_tw,DE_low,'type','spearman');
title(['r: ' num2str(rho_sp,2) ' p=' num2str(p_sp,2)]);
text(-0.19,60,'D','Fontsize',20)


%% check RFT distractor noisy load correlation for mention in paper

t1=dsearchn(cl_time',0.5)+cl_t1-1;
t2=dsearchn(cl_time',1.35)+cl_t1-1;
RFT_LMI_noisy=mean(nanmean(timecourses.RFT.LMI_dist_noisy(:,:,t1:t2),3),2);
DE_low=RT_all(:,3)-RT_all(:,1);
DE_high=RT_all(:,4)-RT_all(:,2);
[rho,p]=corr(RFT_LMI_noisy,DE_low,'type','spearman');
disp(['Noisy distractor RFT correlation with behavioral distractor interference with low load: p=' num2str(p,2)])
[rho,p]=corr(RFT_LMI_noisy,DE_high,'type','spearman');
disp(['Noisy distractor RFT correlation with behavioral distractor interference with high load: p=' num2str(p,2)])

%% Figure 5 - RFT load effects
f5=figure('Name','Fig 5: RFT load effects');

%time window plots
cl_tw=[0.5 1.35];

cl_tw_ind=[dsearchn(0.35+timecourses.time',cl_tw(1)) dsearchn(0.35+timecourses.time',cl_tw(2))];

for i=1:4
    tmp=[timecourses.RFT.dist{:,i}]';
    dist_cl(:,i)=nanmean(tmp(:,cl_tw_ind(1):cl_tw_ind(2)),2);
    
    tmp=[timecourses.RFT.target{:,i}]';
    target_cl(:,i)=nanmean(tmp(:,cl_tw_ind(1):cl_tw_ind(2)),2);
end

%RFT target ANOVA
subplot(4,4,[1 2])
beh_plot=[target_cl(:,1) target_cl(:,3) target_cl(:,2) target_cl(:,4)];

h=violinplot(-100+beh_plot*100,{'Dist. Noisy','Dist. Salient','Dist. Noisy','Dist. Salient'},'ShowMean',true);
h(1,1).ViolinColor=[0.2 0.2 0.8];
h(1,1).EdgeColor=[0.2 0.2 0.8];
h(1,1).ViolinPlot.LineStyle='--';
h(1,1).ViolinPlot.LineWidth=2;
h(1,2).ViolinColor=[0.2 0.2 0.8];
h(1,2).EdgeColor=[0.2 0.2 0.8];
h(1,2).ViolinPlot.LineWidth=2;
h(1,3).ViolinColor=[0.8 0.2 0.2];
h(1,3).EdgeColor=[0.8 0.2 0.2];
h(1,3).ViolinPlot.LineWidth=2;
h(1,3).ViolinPlot.LineStyle='--';
h(1,4).ViolinColor=[0.8 0.2 0.2];
h(1,4).EdgeColor=[0.8 0.2 0.2];
h(1,4).ViolinPlot.LineWidth=2;
xlabel('Low target load                 high target load')
ylabel('RFT power (rel. baseline %change)');ylim([-100 120]);hline(1,'k');
title('Target RFT')

%add ANOVA significance
line([2 4],[-70 -70],'color','k','linewidth',2);
text(3,-85,'*','Fontsize',16);
text(-1,120,'A','Fontsize',20);

%RFT distractor ANOVA
subplot(4,4,[3 4])
beh_plot=[dist_cl(:,1) dist_cl(:,3) dist_cl(:,2) dist_cl(:,4)];

h=violinplot(-100+beh_plot*100,{'Dist. Noisy','Dist. Salient','Dist. Noisy','Dist. Salient'},'ShowMean',true);
h(1,1).ViolinColor=[0.2 0.2 0.8];
h(1,1).EdgeColor=[0.2 0.2 0.8];
h(1,1).ViolinPlot.LineStyle='--';
h(1,1).ViolinPlot.LineWidth=2;
h(1,2).ViolinColor=[0.2 0.2 0.8];
h(1,2).EdgeColor=[0.2 0.2 0.8];
h(1,2).ViolinPlot.LineWidth=2;
h(1,3).ViolinColor=[0.8 0.2 0.2];
h(1,3).EdgeColor=[0.8 0.2 0.2];
h(1,3).ViolinPlot.LineWidth=2;
h(1,3).ViolinPlot.LineStyle='--';
h(1,4).ViolinColor=[0.8 0.2 0.2];
h(1,4).EdgeColor=[0.8 0.2 0.2];
h(1,4).ViolinPlot.LineWidth=2;
xlabel('Low target load                 high target load');
ylim([-100 120]);hline(1,'k');ax=gca;ax.YTickLabel={};
title('Distractor RFT');

%add ANOVA significance
line([1 3],[70 70],'color','k','linewidth',2); %low vs high load noisy dist
text(2,80,'*','Fontsize',16);
line([3 4],[85 85],'color','k','linewidth',2); %salient vs noisy dist high load
text(3.4,95,'*','Fontsize',16);

%Target RFT timecourse
RFT_target=timecourses.RFT.target;
RFT_dt_target=timecourses.RFT.target_dt;

%truncate to non-nan (aka numbers)
cl_t1=find(~isnan(RFT_target{1,1}),1,'first');
cl_t2=find(~isnan(RFT_target{1,1}),1,'last');

cl_time=timecourses.time(cl_t1:cl_t2)+.35; %cue onset = 0

dt_t1=find(~isnan(RFT_dt_target{1,1}),1,'first');
dt_t2=find(~isnan(RFT_dt_target{1,1}),1,'last');

dt_time=timecourses.dt_time(dt_t1:dt_t2);

%get mean and SE per load condition
for l=1:4
    m_RFT_target{l}=mean([RFT_target{:,l}],2);m_RFT_target{l}=m_RFT_target{l}(cl_t1:cl_t2)*100;
    se_RFT_target{l}=std([RFT_target{:,l}]')./sqrt(Nsubs);se_RFT_target{l}=se_RFT_target{l}(cl_t1:cl_t2)*100;
    
    %DT
    m_RFT_dt_target{l}=mean([RFT_dt_target{:,l}],2);m_RFT_dt_target{l}=m_RFT_dt_target{l}(dt_t1:dt_t2)*100;
    se_RFT_dt_target{l}=std([RFT_dt_target{:,l}]')./sqrt(Nsubs);se_RFT_dt_target{l}=se_RFT_dt_target{l}(dt_t1:dt_t2)*100;
end

%Cue-locked
subplot(4,4,[5 6])
h1=fill([cl_time fliplr(cl_time)],[m_RFT_target{3}'-se_RFT_target{3} fliplr(m_RFT_target{3}'+se_RFT_target{3})],[0 0 1],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([cl_time fliplr(cl_time)],[m_RFT_target{4}'-se_RFT_target{4} fliplr(m_RFT_target{4}'+se_RFT_target{4})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(cl_time,m_RFT_target{3},'color',[0 0 1],'linewidth',2);
h4=plot(cl_time,m_RFT_target{4},'color',[1 0 0],'linewidth',2);
xlim([-0.25 cl_time(end)]);ylim([55 110]);vline(0,'color','k','linestyle','--','linewidth',2);
if timecourse_stats, plot_stat(RFT_target_cl_salient_stat,0.6,[0 0 0]); end
text(0.05,60,'CUE onset');title('Target RFT (salient distractors)');legend([h3,h4],{'low target load','high target load'},'location','southeast')
xlabel('Time(s)');ylabel('RFT power (%Baseline)')
text(-0.5,110,'B','Fontsize',20)

%discrimination target-locked
subplot(4,4,[7 8])
h1=fill([dt_time fliplr(dt_time)],[m_RFT_dt_target{3}'-se_RFT_dt_target{3} fliplr(m_RFT_dt_target{3}'+se_RFT_dt_target{3})],[0 0 1],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([dt_time fliplr(dt_time)],[m_RFT_dt_target{4}'-se_RFT_dt_target{4} fliplr(m_RFT_dt_target{4}'+se_RFT_dt_target{4})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(dt_time,m_RFT_dt_target{3},'color',[0 0 1],'linewidth',2);
h4=plot(dt_time,m_RFT_dt_target{4},'color',[1 0 0],'linewidth',2);
xlim([-cl_time(end) 0.25]);ylim([55 110]);vline(0,'color','k','linestyle','--','linewidth',2);
if timecourse_stats, plot_stat(RFT_target_dt_salient_stat,0.6,[0 0 0]); end
text(-0.55,60,'Discrimination target');
xlabel('Time(s)');ax=gca;ax.YTickLabel={};

%Add empty space for source loc
subplot(4,4,[9 10])
text(-0.15,1,'C','Fontsize',20)
axis off;

%correlations
DE_time=timecourses.time+0.35;
t1=dsearchn(DE_time',0.5);
t2=dsearchn(DE_time',1.35);
RFT_DE_low_tw=nanmean(nanmean(timecourses.RFT.DE_low_load(:,:,t1:t2),2),3);
RFT_DE_high_tw=nanmean(nanmean(timecourses.RFT.DE_high_load(:,:,t1:t2),2),3);
DE_high=RT_all(:,4)-RT_all(:,2);

%Low load
subplot(4,4,11)
scatter(RFT_DE_low_tw,DE_low);axis([-0.08 0.08 -50 60])
h=lsline;set(h,'color','r');hline(0,'k');vline(0,'k');
xlabel('Target RFT - Distractor effect');ylabel('Beh. dist. interference (ms)')
[rho_sp,p_sp]=corr(RFT_DE_low_tw,DE_low,'type','spearman');
title({'Low target load',['r: ' num2str(rho_sp,3) ' p=' num2str(p_sp,2)]})
text(-0.135,60,'D','Fontsize',20)

%High load
subplot(4,4,12)
scatter(RFT_DE_high_tw,DE_high);axis([-0.08 0.08 -50 60])
h=lsline;set(h,'color','r');hline(0,'k');vline(0,'k');
xlabel('Target RFT - Distractor effect');
[rho_sp,p_sp]=corr(RFT_DE_high_tw,DE_high,'type','spearman');
title({'High target load',['r: ' num2str(rho_sp,3) ' p=' num2str(p_sp,2)]})
ax=gca;ax.YTickLabel={};

%Distractor RFT timecourse noisy low vs high
RFT_dist=timecourses.RFT.dist;
RFT_dt_dist=timecourses.RFT.dist_dt;

%Get MEANs and SEs of the timecourses
Nsubs=size(RFT_dist,1);

%get mean and SE per load condition
for l=1:4
    
    %RFT
    m_RFT_distractor{l}=mean([RFT_dist{selection,l}],2);m_RFT_distractor{l}=m_RFT_distractor{l}(cl_t1:cl_t2)*100;
    se_RFT_distractor{l}=std([RFT_dist{selection,l}]')./sqrt(Nsubs);se_RFT_distractor{l}=se_RFT_distractor{l}(cl_t1:cl_t2)*100;
    
    %RFT DT
    m_RFT_dt_distractor{l}=mean([RFT_dt_dist{selection,l}],2);m_RFT_dt_distractor{l}=m_RFT_dt_distractor{l}(dt_t1:dt_t2)*100;
    se_RFT_dt_distractor{l}=std([RFT_dt_dist{selection,l}]')./sqrt(Nsubs);se_RFT_dt_distractor{l}=se_RFT_dt_distractor{l}(dt_t1:dt_t2)*100;
end

%Cue-locked
subplot(4,4,[13 14])
h1=fill([cl_time fliplr(cl_time)],[m_RFT_distractor{1}'-se_RFT_distractor{1} fliplr(m_RFT_distractor{1}'+se_RFT_distractor{1})],[0 0 1],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([cl_time fliplr(cl_time)],[m_RFT_distractor{2}'-se_RFT_distractor{2} fliplr(m_RFT_distractor{2}'+se_RFT_distractor{2})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(cl_time,m_RFT_distractor{1},'color',[0 0 1],'linewidth',2,'linestyle','--');
h4=plot(cl_time,m_RFT_distractor{2},'color',[1 0 0],'linewidth',2,'linestyle','--');
xlim([-0.25 cl_time(end)]);ylim([45 95]);vline(0,'color','k','linestyle','--','linewidth',2);
text(0.05,48,'CUE onset');title('Distractor RFT (noisy distractors)')
xlabel('Time(s)');ylabel('RFT power (%Baseline)')
legend([h3,h4],{'low target load','high target load'},'location','northeast')
text(-0.45,95,'E','Fontsize',20)

%Discrimination target-locked
subplot(4,4,[15 16])
h1=fill([dt_time fliplr(dt_time)],[m_RFT_dt_distractor{1}'-se_RFT_dt_distractor{1} fliplr(m_RFT_dt_distractor{1}'+se_RFT_dt_distractor{1})],[0 0 1],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([dt_time fliplr(dt_time)],[m_RFT_dt_distractor{2}'-se_RFT_dt_distractor{2} fliplr(m_RFT_dt_distractor{2}'+se_RFT_dt_distractor{2})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(dt_time,m_RFT_dt_distractor{1},'color',[0 0 1],'linewidth',2,'linestyle','--');
h4=plot(dt_time,m_RFT_dt_distractor{2},'color',[1 0 0],'linewidth',2,'linestyle','--');
xlim([-cl_time(end) 0.25]);ylim([45 95]);vline(0,'color','k','linestyle','--','linewidth',2);
text(-0.45,48,'Discrimination target');
xlabel('Time(s)');ax=gca;ax.YTickLabel={};

%% Figure 6 - Effects of distractor salience
f6=figure('Name','Fig 6: Salience effects');

%alpha
cl_t1=find(~isnan(Alpha_dist{1,1}),1,'first');
cl_t2=find(~isnan(Alpha_dist{1,1}),1,'last');

cl_time=timecourses.time(cl_t1:cl_t2)+.35; %cue onset = 0

dt_t1=find(~isnan(Alpha_dt_dist{1,1}),1,'first');
dt_t2=find(~isnan(Alpha_dt_dist{1,1}),1,'last');

dt_time=timecourses.dt_time(dt_t1:dt_t2);

%Cue-locked Distractor Alpha salience
subplot(2,4,[1 2])
h1=fill([cl_time fliplr(cl_time)],[m_Alpha_distractor{2}'-se_Alpha_distractor{2} fliplr(m_Alpha_distractor{2}'+se_Alpha_distractor{2})],[1 0 0],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([cl_time fliplr(cl_time)],[m_Alpha_distractor{4}'-se_Alpha_distractor{4} fliplr(m_Alpha_distractor{4}'+se_Alpha_distractor{4})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(cl_time,m_Alpha_distractor{2},'color',[1 0 0],'linewidth',2,'linestyle','--');
h4=plot(cl_time,m_Alpha_distractor{4},'color',[1 0 0],'linewidth',2);
xlim([-0.25 cl_time(end)]);ylim([93 138]);vline(0,'color','k','linestyle','--','linewidth',2);
if timecourse_stats, plot_stat(Alpha_dist_cl_highload_stat,0.97,[0 0 0]); end
text(0.05,96,'CUE onset');title('Distractor Alpha (high load)')
xlabel('Time(s)');ylabel('Rel. Baseline Alpha power')
legend([h3,h4],{'Noisy distractor','Salient distractor'},'location','northeast')
text(-0.6,138,'A','Fontsize',20)

%Discrimination target-locked
subplot(2,4,[3 4])
h1=fill([dt_time fliplr(dt_time)],[m_Alpha_dt_distractor{2}'-se_Alpha_dt_distractor{2} fliplr(m_Alpha_dt_distractor{2}'+se_Alpha_dt_distractor{1})],[1 0 0],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([dt_time fliplr(dt_time)],[m_Alpha_dt_distractor{4}'-se_Alpha_dt_distractor{4} fliplr(m_Alpha_dt_distractor{4}'+se_Alpha_dt_distractor{3})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(dt_time,m_Alpha_dt_distractor{2},'color',[1 0 0],'linewidth',2,'linestyle','--');
h4=plot(dt_time,m_Alpha_dt_distractor{4},'color',[1 0 0],'linewidth',2);
xlim([-cl_time(end) 0.25]);ylim([93 138]);vline(0,'color','k','linestyle','--','linewidth',2);
if timecourse_stats, plot_stat(Alpha_dist_dt_highload_stat,0.97,[0 0 0]); end
text(-0.6,96,'Discrimination target');
xlabel('Time(s)');ax=gca;ax.YTickLabel={};

%RFT Distractor salience
%truncate to non-nan (aka numbers)
cl_t1=find(~isnan(RFT_target{1,1}),1,'first');
cl_t2=find(~isnan(RFT_target{1,1}),1,'last');

cl_time=timecourses.time(cl_t1:cl_t2)+.35; %cue onset = 0

dt_t1=find(~isnan(RFT_dt_target{1,1}),1,'first');
dt_t2=find(~isnan(RFT_dt_target{1,1}),1,'last');

dt_time=timecourses.dt_time(dt_t1:dt_t2);

%Cue-locked
subplot(2,4,[5 6])
h1=fill([cl_time fliplr(cl_time)],[m_RFT_distractor{2}'-se_RFT_distractor{2} fliplr(m_RFT_distractor{2}'+se_RFT_distractor{2})],[1 0 0],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([cl_time fliplr(cl_time)],[m_RFT_distractor{4}'-se_RFT_distractor{4} fliplr(m_RFT_distractor{4}'+se_RFT_distractor{4})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(cl_time,m_RFT_distractor{2},'color',[1 0 0],'linewidth',2,'linestyle','--');
h4=plot(cl_time,m_RFT_distractor{4},'color',[1 0 0],'linewidth',2);
xlim([-0.25 cl_time(end)]);ylim([45 95]);vline(0,'color','k','linestyle','--','linewidth',2);
if timecourse_stats, plot_stat(RFT_dist_highload_stat,0.47,[0 0 0]); end
text(0.05,48,'CUE onset');title('Distractor RFT (high load)')
xlabel('Time(s)');ylabel('Rel. Baseline RFT power')
text(-0.6,95,'B','Fontsize',20)

%Discrimination target-locked
subplot(2,4,[7 8])
h1=fill([dt_time fliplr(dt_time)],[m_RFT_dt_distractor{2}'-se_RFT_dt_distractor{2} fliplr(m_RFT_dt_distractor{2}'+se_RFT_dt_distractor{2})],[1 0 0],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([dt_time fliplr(dt_time)],[m_RFT_dt_distractor{4}'-se_RFT_dt_distractor{4} fliplr(m_RFT_dt_distractor{4}'+se_RFT_dt_distractor{4})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(dt_time,m_RFT_dt_distractor{2},'color',[1 0 0],'linewidth',2,'linestyle','--');
h4=plot(dt_time,m_RFT_dt_distractor{4},'color',[1 0 0],'linewidth',2);
xlim([-cl_time(end) 0.25]);ylim([45 95]);vline(0,'color','k','linestyle','--','linewidth',2);
if timecourse_stats, plot_stat(RFT_dist_dt_highload_stat,0.47,[0 0 0]); end
text(-0.6,48,'Discrimination target');legend([h3,h4],{'Noisy distractors','Salient distractors'},'location','northwest')
xlabel('Time(s)');ax=gca;ax.YTickLabel={};

%% Calculate correlations with behaviour for paper (mention)
DE_dist_high=(([timecourses.alpha.dist{:,4}]-[timecourses.alpha.dist{:,2}])./([timecourses.alpha.dist{:,4}]+[timecourses.alpha.dist{:,2}]))';
DE_time=timecourses.time+0.35;
t1=dsearchn(DE_time',0.5);
t2=dsearchn(DE_time',1.35);

%Alpha
Alpha_DE_dist_high_tw=nanmean(DE_dist_high(:,t1:t2),2);
[rho,p]=corr(Alpha_DE_dist_high_tw,DE_high);
disp(['High load Distractor alpha salient vs noisy correlation with behavioral distractor interference with highload: r: ' num2str(rho,2) ' p=' num2str(p,2) ' (pearson)'])
[rho,p]=corr(Alpha_DE_dist_high_tw,DE_high,'type','spearman');
disp(['High load Distractor alpha salient vs noisy correlation with behavioral distractor interference with highload: r: ' num2str(rho,2) ' p=' num2str(p,2) ' (spearman)'])

%RFT
RFT_DE_dist_high=(([timecourses.RFT.dist{:,4}]-[timecourses.RFT.dist{:,2}])./([timecourses.RFT.dist{:,4}]+[timecourses.RFT.dist{:,2}]))';
RFT_DE_dist_high_tw=nanmean(RFT_DE_dist_high(:,t1:t2),2);
[rho,p]=corr(RFT_DE_dist_high_tw,DE_high);
disp(['High load Distractor RFT salient vs noisy correlation with behavioral distractor interference with highload: r: ' num2str(rho,2) ' p=' num2str(p,2) ' (pearson)'])
[rho,p]=corr(RFT_DE_dist_high_tw,DE_high,'type','spearman');
disp(['High load Distractor RFT salient vs noisy correlation with behavioral distractor interference with highload: r: ' num2str(rho,2) ' p=' num2str(p,2) ' (spearman)'])



