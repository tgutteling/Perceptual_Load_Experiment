function [] = D2_Paper_figures()
%
% Here we only plot all the relevant figures as reported in the paper
%

%folders
group_folder='Z:\Tjerk\Load2\proc_data\group\';

%https://github.com/bastibe/Violinplot-Matlab/blob/master/Violin.m
addpath Z:\Tjerk\Scripts_general\Violinplot\

%vline and hline
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
h=bar([mean(beh_plot(:,1)) mean(beh_plot(:,2)) ; mean(beh_plot(:,3)) mean(beh_plot(:,4))],'FaceColor','flat','FaceAlpha',0.5,'EdgeColor','flat','LineWidth',1.5);
hold on;
errorbar([h(1).XEndPoints h(2).XEndPoints],[h(1).YData h(2).YData],std(beh_plot)/sqrt(length(beh_plot)),'.','Color','k','LineWidth',1.5)
ylim([250 530])
h(1).CData(1,:)=[0.2 0.2 0.8];
h(1).CData(2,:)=[0.8 0.2 0.2];
h(2).CData(1,:)=[0.2 0.2 0.8];
h(2).CData(2,:)=[0.8 0.2 0.2];
h(1).LineStyle='--';
ylabel('Reaction time (ms)');
%xticklabels({'Low target load','High target load'})
xticklabels({'noisy distractors   salient distractors','noisy distractors   salient distractors'})
xlabel('Low target load                               High target load')
text(0.2,530,'A','Fontsize',20);box off

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
subplot(5,6,1:3)
beh_plot=[target_cl(:,3) target_cl(:,4) target_cl(:,1) target_cl(:,2)]*100;

h=bar([mean(beh_plot(:,1)) mean(beh_plot(:,2)) ; mean(beh_plot(:,3)) mean(beh_plot(:,4))],'FaceColor','flat','FaceAlpha',0.5,'EdgeColor','flat','LineWidth',1.5);
hold on;
errorbar([h(1).XEndPoints h(2).XEndPoints],[h(1).YData h(2).YData],std(beh_plot)/sqrt(length(beh_plot)),'.','Color','k','LineWidth',1.5)
ylim([90 130])
h(1).CData(1,:)=[0.2 0.2 0.8];
h(1).CData(2,:)=[0.2 0.2 0.8];
h(2).CData(1,:)=[0.8 0.2 0.2];
h(2).CData(2,:)=[0.8 0.2 0.2];
%h(1).LineStyle='--';
%ylabel('Alpha power (rel. baseline, %change)');hline(0,'k');
ylabel('Alpha power (%baseline)');hline(0,'k');
xticklabels({'Low load  high load','low load  high load'})
xlabel('Salient distractors     Noisy distractors')
title('Target Alpha')
text(0.1,130,'A','Fontsize',20);box off

%Alpha ANOVA results DISTRACTOR
subplot(5,6,4:6)
beh_plot=[dist_cl(:,3) dist_cl(:,4) dist_cl(:,1) dist_cl(:,2)]*100;        
     
h=bar([mean(beh_plot(:,1)) mean(beh_plot(:,2)) ; mean(beh_plot(:,3)) mean(beh_plot(:,4))],'FaceColor','flat','FaceAlpha',0.5,'EdgeColor','flat','LineWidth',1.5);
hold on;
errorbar([h(1).XEndPoints h(2).XEndPoints],[h(1).YData h(2).YData],std(beh_plot)/sqrt(length(beh_plot)),'.','Color','k','LineWidth',1.5)
ylim([90 130])
h(1).CData(1,:)=[0.2 0.2 0.8];
h(1).CData(2,:)=[0.2 0.2 0.8];
h(2).CData(1,:)=[0.8 0.2 0.2];
h(2).CData(2,:)=[0.8 0.2 0.2];
%h(1).LineStyle='--';
xticklabels({'Low load  high load','low load  high load'})
xlabel('Salient distractors     Noisy distractors')
title('Distractor Alpha');box off

%add ANOVA significance
line([h(1).XEndPoints(1) h(2).XEndPoints(1)],[121 121],'color','k','linewidth',2) %low vs high load salient dist
line([h(1).XEndPoints(1) h(1).XEndPoints(1)],[115 121],'color','k','linewidth',2)
line([h(2).XEndPoints(1) h(2).XEndPoints(1)],[120 121],'color','k','linewidth',2)
line(repmat(mean([h(1).XEndPoints(1) h(2).XEndPoints(1)]),1,2),[121 122],'color','k','linewidth',2)
text(mean([h(1).XEndPoints(1) h(2).XEndPoints(1)]),122.5,'*','Fontsize',16,'HorizontalAlignment','center')

line([h(2).XEndPoints(1) h(2).XEndPoints(2)],[123 123],'color','k','linewidth',2) %high load noisy vs salient dist
line([h(2).XEndPoints(1) h(2).XEndPoints(1)],[122 123],'color','k','linewidth',2)
line([h(2).XEndPoints(2) h(2).XEndPoints(2)],[115 123],'color','k','linewidth',2)
line(repmat(mean([h(2).XEndPoints(1) h(2).XEndPoints(2)]),1,2),[123 124],'color','k','linewidth',2)
text(mean([h(2).XEndPoints(1)  h(2).XEndPoints(2)]),124.5,'*','Fontsize',16,'HorizontalAlignment','center')

%SECTION B - Zoom in on significant effects
%Load effect (distractor alpha)
subplot(5,6,7:8)
beh_plot=[dist_cl(:,4)-dist_cl(:,3) dist_cl(:,2)-dist_cl(:,1)]*100;        
h=violinplot(beh_plot,{'Salient Distractors','Noisy Distractors'},'ShowMean',true);
h(1,1).ViolinColor=[0 0 0];
h(1,1).EdgeColor=[0 0 0];
h(1,1).ViolinPlot.LineWidth=2;
h(1,2).ViolinColor=[0 0 0];
h(1,2).EdgeColor=[0 0 0];
h(1,2).ViolinPlot.LineStyle='--';
h(1,2).ViolinPlot.LineWidth=2;
ylabel('\DeltaDistractor Alpha')
title('Distractor Alpha - load effect')
text(1,40,'*','Fontsize',30,'HorizontalAlignment','center')
ylim([-50 50])
text(0.01,50,'B','Fontsize',20);hline(0,'k');

%SECTION C - Salience effect (distractor alpha)
subplot(5,6,10:11)
beh_plot=[dist_cl(:,4)-dist_cl(:,2) dist_cl(:,3)-dist_cl(:,1)]*100;        
h=violinplot(beh_plot,{'High load','Low load'},'ShowMean',true);
h(1,1).ViolinColor=[0.8 0.2 0.2];
h(1,1).EdgeColor=[0.8 0.2 0.2];
h(1,1).ViolinPlot.LineWidth=2;
h(1,2).ViolinColor=[0.2 0.2 0.8];
h(1,2).EdgeColor=[0.2 0.2 0.8];
h(1,2).ViolinPlot.LineWidth=2;
ylabel('\DeltaDistractor Alpha')
title('Distractor Alpha - salience effect')
text(1,40,'*','Fontsize',30,'HorizontalAlignment','center')
ylim([-50 50]);hline(0,'k')
text(0.01,50,'C','Fontsize',20);hline(0,'k');


%SECTION D - Alpha distractor time course
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

subplot(5,6,13:15)
h1=fill([cl_time fliplr(cl_time)],[m_Alpha_distractor{3}'-se_Alpha_distractor{3} fliplr(m_Alpha_distractor{3}'+se_Alpha_distractor{3})],[0 0 1],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([cl_time fliplr(cl_time)],[m_Alpha_distractor{4}'-se_Alpha_distractor{4} fliplr(m_Alpha_distractor{4}'+se_Alpha_distractor{4})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(cl_time,m_Alpha_distractor{3},'color',[0 0 1],'linewidth',2);
h4=plot(cl_time,m_Alpha_distractor{4},'color',[1 0 0],'linewidth',2);
xlim([-0.25 cl_time(end)]);ylim([90 140]);vline(0,'color','k','linestyle','--','linewidth',2);
text(0.05,94,'CUE onset');title('Distractor Alpha (salient distractors)')
xlabel('Time(s)');ylabel('Alpha power (%Baseline)')
legend([h3,h4],{'low target load','high target load'},'location','northeast');box off
text(-0.57,140,'D','Fontsize',20)

subplot(5,6,16:18)
h1=fill([dt_time fliplr(dt_time)],[m_Alpha_dt_distractor{3}'-se_Alpha_dt_distractor{3} fliplr(m_Alpha_dt_distractor{3}'+se_Alpha_dt_distractor{3})],[0 0 1],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([dt_time fliplr(dt_time)],[m_Alpha_dt_distractor{4}'-se_Alpha_dt_distractor{4} fliplr(m_Alpha_dt_distractor{4}'+se_Alpha_dt_distractor{4})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(dt_time,m_Alpha_dt_distractor{3},'color',[0 0 1],'linewidth',2);
h4=plot(dt_time,m_Alpha_dt_distractor{4},'color',[1 0 0],'linewidth',2);
xlim([-cl_time(end) 0.25]);ylim([90 140]);vline(0,'color','k','linestyle','--','linewidth',2);
text(-0.55,94,'Discrimination target');
xlabel('Time(s)');ax=gca;ax.YTickLabel={};box off

%SECTION E- source loc (added externally)
subplot(5,6,19:21)
text(-0.15,1,'E','Fontsize',20)
box off;axis off

%SECTION F - correlation
subplot(5,6,22:23)
text(-0.15,1,'F','Fontsize',20)
t1=dsearchn(cl_time',0.5)+cl_t1-1;
t2=dsearchn(cl_time',1.35)+cl_t1-1;
Alpha_LMI_salient_tw=mean(timecourses.alpha.LMI_dist_salient(:,t1:t2),2);
DE_low=RT_all(:,3)-RT_all(:,1);
scatter(Alpha_LMI_salient_tw,DE_low);
h=lsline;set(h,'color','r');hline(0,'k');vline(0,'k');
xlabel('Distractor Alpha - Load modulation index');ylabel('Beh. low load dist. interference (ms)')
[rho_sp,p_sp]=corr(Alpha_LMI_salient_tw,DE_low,'type','spearman');
title(['r: ' num2str(rho_sp,2) ' p=' num2str(p_sp,2)]); 
text(-0.19,60,'F','Fontsize',20)

%SECTION G - Distractor salience timecourse
subplot(5,6,[25:27])
h1=fill([cl_time fliplr(cl_time)],[m_Alpha_distractor{2}'-se_Alpha_distractor{2} fliplr(m_Alpha_distractor{2}'+se_Alpha_distractor{2})],[1 0 0],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([cl_time fliplr(cl_time)],[m_Alpha_distractor{4}'-se_Alpha_distractor{4} fliplr(m_Alpha_distractor{4}'+se_Alpha_distractor{4})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(cl_time,m_Alpha_distractor{2},'color',[1 0 0],'linewidth',2,'linestyle','--');
h4=plot(cl_time,m_Alpha_distractor{4},'color',[1 0 0],'linewidth',2);
xlim([-0.25 cl_time(end)]);ylim([93 138]);vline(0,'color','k','linestyle','--','linewidth',2);
text(0.05,96,'CUE onset');title('Distractor Alpha (high load)')
xlabel('Time(s)');ylabel('Alpha power (%Baseline)')
legend([h3,h4],{'Noisy distractor','Salient distractor'},'location','northeast'); box off
text(-0.6,138,'G','Fontsize',20)
subplot(5,6,[28:30])
h1=fill([dt_time fliplr(dt_time)],[m_Alpha_dt_distractor{2}'-se_Alpha_dt_distractor{2} fliplr(m_Alpha_dt_distractor{2}'+se_Alpha_dt_distractor{1})],[1 0 0],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([dt_time fliplr(dt_time)],[m_Alpha_dt_distractor{4}'-se_Alpha_dt_distractor{4} fliplr(m_Alpha_dt_distractor{4}'+se_Alpha_dt_distractor{3})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(dt_time,m_Alpha_dt_distractor{2},'color',[1 0 0],'linewidth',2,'linestyle','--');
h4=plot(dt_time,m_Alpha_dt_distractor{4},'color',[1 0 0],'linewidth',2);
xlim([-cl_time(end) 0.25]);ylim([93 138]);vline(0,'color','k','linestyle','--','linewidth',2);
text(-0.6,96,'Discrimination target');
xlabel('Time(s)');ax=gca;ax.YTickLabel={};box off


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

%% Figure 5 - target RFT 
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

%SECTION A - timewindow ANOVA 
subplot(3,8,1:4)
beh_plot=[target_cl(:,3) target_cl(:,4) target_cl(:,1) target_cl(:,2)]*100;

h=bar([mean(beh_plot(:,1)) mean(beh_plot(:,2)) ; mean(beh_plot(:,3)) mean(beh_plot(:,4))],'FaceColor','flat','FaceAlpha',0.5,'EdgeColor','flat','LineWidth',1.5);
hold on;
errorbar([h(1).XEndPoints h(2).XEndPoints],[h(1).YData h(2).YData],std(beh_plot)/sqrt(length(beh_plot)),'.','Color','k','LineWidth',1.5)
ylim([65 100])
h(1).CData(1,:)=[0.2 0.2 0.8];
h(1).CData(2,:)=[0.2 0.2 0.8];
h(2).CData(1,:)=[0.8 0.2 0.2];
h(2).CData(2,:)=[0.8 0.2 0.2];
ylabel('RFT power (%baseline)');
xticklabels({'Low load  high load','low load  high load'})
xlabel('Salient distractors     Noisy distractors')
title('Target RFT')
text(0.1,100,'A','Fontsize',20);box off

%add ANOVA significance
line([h(1).XEndPoints(1) h(2).XEndPoints(1)],[95 95],'color','k','linewidth',2) %low vs high load salient dist
line([h(1).XEndPoints(1) h(1).XEndPoints(1)],[95 85],'color','k','linewidth',2)
line([h(2).XEndPoints(1) h(2).XEndPoints(1)],[95 93],'color','k','linewidth',2)
line(repmat(mean([h(1).XEndPoints(1) h(2).XEndPoints(1)]),1,2),[95 97],'color','k','linewidth',2)
text(mean([h(1).XEndPoints(1) h(2).XEndPoints(1)]),98,'*','Fontsize',16,'HorizontalAlignment','center')

%SECTION B - ZOOM IN ON LOAD EFFECT 
subplot(3,8,6:8)
beh_plot=[target_cl(:,4)-target_cl(:,3) target_cl(:,2)-target_cl(:,1)]*100;        
h=violinplot(beh_plot,{'Salient Distractors','Noisy Distractors'},'ShowMean',true);
h(1,1).ViolinColor=[0 0 0];
h(1,1).EdgeColor=[0 0 0];
h(1,1).ViolinPlot.LineWidth=2;
h(1,2).ViolinColor=[0 0 0];
h(1,2).EdgeColor=[0 0 0];
h(1,2).ViolinPlot.LineStyle='--';
h(1,2).ViolinPlot.LineWidth=2;
ylabel('\DeltaTarget RFT')
title('Target load effect')
text(1,30,'*','Fontsize',30,'HorizontalAlignment','center')
ylim([-25 40]);
text(0.01,40,'B','Fontsize',20);hline(0,'k');

%SECTION C - TARGET RFT TIMECOURSE 
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

subplot(3,8,9:12)
h1=fill([cl_time fliplr(cl_time)],[m_RFT_target{3}'-se_RFT_target{3} fliplr(m_RFT_target{3}'+se_RFT_target{3})],[0 0 1],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([cl_time fliplr(cl_time)],[m_RFT_target{4}'-se_RFT_target{4} fliplr(m_RFT_target{4}'+se_RFT_target{4})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(cl_time,m_RFT_target{3},'color',[0 0 1],'linewidth',2);
h4=plot(cl_time,m_RFT_target{4},'color',[1 0 0],'linewidth',2);
xlim([-0.25 cl_time(end)]);ylim([55 110]);vline(0,'color','k','linestyle','--','linewidth',2);
text(0.05,60,'CUE onset');title('Target RFT (salient distractors)');legend([h3,h4],{'low target load','high target load'},'location','southeast')
xlabel('Time(s)');ylabel('RFT power (%Baseline)');box off
text(-0.5,110,'C','Fontsize',20)

subplot(3,8,13:16)
h1=fill([dt_time fliplr(dt_time)],[m_RFT_dt_target{3}'-se_RFT_dt_target{3} fliplr(m_RFT_dt_target{3}'+se_RFT_dt_target{3})],[0 0 1],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([dt_time fliplr(dt_time)],[m_RFT_dt_target{4}'-se_RFT_dt_target{4} fliplr(m_RFT_dt_target{4}'+se_RFT_dt_target{4})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(dt_time,m_RFT_dt_target{3},'color',[0 0 1],'linewidth',2);
h4=plot(dt_time,m_RFT_dt_target{4},'color',[1 0 0],'linewidth',2);
xlim([-cl_time(end) 0.25]);ylim([55 110]);vline(0,'color','k','linestyle','--','linewidth',2);
text(-0.55,60,'Discrimination target');
xlabel('Time(s)');ax=gca;ax.YTickLabel={};box off

%SECTION D - SOURCE LOCALIZATION
%source loc (empty)
subplot(3,8,17:20)
text(-0.15,1,'D','Fontsize',20)
axis off;

%SECTION E - CORRELATIONS
%correlations
DE_time=timecourses.time+0.35;
t1=dsearchn(DE_time',0.5);
t2=dsearchn(DE_time',1.35);
RFT_DE_low_tw=nanmean(nanmean(timecourses.RFT.DE_low_load(:,:,t1:t2),2),3);
RFT_DE_high_tw=nanmean(nanmean(timecourses.RFT.DE_high_load(:,:,t1:t2),2),3);
DE_high=RT_all(:,4)-RT_all(:,2);

subplot(3,8,21:22)
scatter(RFT_DE_low_tw,DE_low);axis([-0.08 0.08 -50 60])
h=lsline;set(h,'color','r');hline(0,'k');vline(0,'k');
xlabel('Target RFT - Distractor effect');ylabel('Beh. dist. interference (ms)')
[rho_sp,p_sp]=corr(RFT_DE_low_tw,DE_low,'type','spearman');
title({'Low target load',['r: ' num2str(rho_sp,3) ' p=' num2str(p_sp,2)]})
text(-0.135,60,'E','Fontsize',20)
subplot(3,8,23:24)
scatter(RFT_DE_high_tw,DE_high);axis([-0.08 0.08 -50 60])
h=lsline;set(h,'color','r');hline(0,'k');vline(0,'k');
xlabel('Target RFT - Distractor effect');
[rho_sp,p_sp]=corr(RFT_DE_high_tw,DE_high,'type','spearman');
title({'High target load',['r: ' num2str(rho_sp,3) ' p=' num2str(p_sp,2)]})
ax=gca;ax.YTickLabel={};

%% Figure 6 - Distractor RFT

f6=figure('Name','Fig 6: Salience effects');
%SECTION A - timewindow ANOVA 
subplot(4,8,1:4)
beh_plot=[dist_cl(:,3) dist_cl(:,4) dist_cl(:,1) dist_cl(:,2)]*100;

h=bar([mean(beh_plot(:,1)) mean(beh_plot(:,2)) ; mean(beh_plot(:,3)) mean(beh_plot(:,4))],'FaceColor','flat','FaceAlpha',0.5,'EdgeColor','flat','LineWidth',1.5);
hold on;
errorbar([h(1).XEndPoints h(2).XEndPoints],[h(1).YData h(2).YData],std(beh_plot)/sqrt(length(beh_plot)),'.','Color','k','LineWidth',1.5)
ylim([40 90])
h(1).CData(1,:)=[0.2 0.2 0.8];
h(1).CData(2,:)=[0.2 0.2 0.8];
h(2).CData(1,:)=[0.8 0.2 0.2];
h(2).CData(2,:)=[0.8 0.2 0.2];
%h(1).LineStyle='--';
ylabel('RFT power (%baseline)');
xticklabels({'Low load  high load','low load  high load'})
xlabel('Salient distractors     Noisy distractors')
title('Distractor RFT')
text(0.1,90,'A','Fontsize',20);box off

%add ANOVA significance
line([h(1).XEndPoints(2) h(2).XEndPoints(2)],[77 77],'color','k','linewidth',2) %low vs high load noisy dist
line([h(1).XEndPoints(2) h(1).XEndPoints(2)],[70 77],'color','k','linewidth',2)
line([h(2).XEndPoints(2) h(2).XEndPoints(2)],[75 77],'color','k','linewidth',2)
line(repmat(mean([h(1).XEndPoints(2) h(2).XEndPoints(2)]),1,2),[77 79],'color','k','linewidth',2)
text(mean([h(1).XEndPoints(2) h(2).XEndPoints(2)]),80,'*','Fontsize',16,'HorizontalAlignment','center')

line([h(2).XEndPoints(1) h(2).XEndPoints(2)],[83 83],'color','k','linewidth',2) %Noisy vs salient high load
line([h(2).XEndPoints(1) h(2).XEndPoints(1)],[75 83],'color','k','linewidth',2)
line([h(2).XEndPoints(2) h(2).XEndPoints(2)],[81 83],'color','k','linewidth',2)
line(repmat(mean([h(2).XEndPoints(1) h(2).XEndPoints(2)]),1,2),[83 85],'color','k','linewidth',2)
text(mean([h(2).XEndPoints(1) h(2).XEndPoints(2)]),86,'*','Fontsize',16,'HorizontalAlignment','center')

%SECTION B - POST-HOC EFFECTS
%zoom in on load effect
%subplot(4,8,9:12)
subplot(4,8,6:8)
beh_plot=[dist_cl(:,4)-dist_cl(:,3) dist_cl(:,2)-dist_cl(:,1)]*100;        
h=violinplot(beh_plot,{'Salient Distractors','Noisy Distractors'},'ShowMean',true);
h(1,1).ViolinColor=[0 0 0];
h(1,1).EdgeColor=[0 0 0];
h(1,1).ViolinPlot.LineWidth=2;
h(1,2).ViolinColor=[0 0 0];
h(1,2).EdgeColor=[0 0 0];
h(1,2).ViolinPlot.LineStyle='--';
h(1,2).ViolinPlot.LineWidth=2;
ylabel('\DeltaDistractor RFT')
title('Distractor RFT - load effect')
text(2,25,'*','Fontsize',30,'HorizontalAlignment','center')
ylim([-40 40]);
text(0.1,40,'B','Fontsize',20);hline(0,'k');

%SECTION C - zoom in on salience effect
subplot(4,8,9:11)
beh_plot=[dist_cl(:,4)-dist_cl(:,2) dist_cl(:,3)-dist_cl(:,1)]*100;        
h=violinplot(beh_plot,{'High Load','Low Load'},'ShowMean',true);
h(1,1).ViolinColor=[0.8 0.2 0.2];
h(1,1).EdgeColor=[0.8 0.2 0.2];
h(1,1).ViolinPlot.LineWidth=2;
h(1,2).ViolinColor=[0.2 0.2 0.8];
h(1,2).EdgeColor=[0.2 0.2 0.8];
h(1,2).ViolinPlot.LineWidth=2;
title('Distractor RFT - salience effect')
text(1,20,'*','Fontsize',30,'HorizontalAlignment','center')
ylabel('\DeltaDistractor RFT')
ylim([-40 40]);hline(0,'k');
text(0.1,40,'C','Fontsize',20);

%SECTION D - SOURCE LOC (EMPTY)
subplot(4,8,12:16)
text(0,1,'D','Fontsize',20)
axis off;

%SECTION E - distractor RFT load effect timecourse
%Distractor RFT noisy low vs high
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

subplot(4,8,17:20)
h1=fill([cl_time fliplr(cl_time)],[m_RFT_distractor{1}'-se_RFT_distractor{1} fliplr(m_RFT_distractor{1}'+se_RFT_distractor{1})],[0 0 1],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([cl_time fliplr(cl_time)],[m_RFT_distractor{2}'-se_RFT_distractor{2} fliplr(m_RFT_distractor{2}'+se_RFT_distractor{2})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(cl_time,m_RFT_distractor{1},'color',[0 0 1],'linewidth',2,'linestyle','--');
h4=plot(cl_time,m_RFT_distractor{2},'color',[1 0 0],'linewidth',2,'linestyle','--');
xlim([-0.25 cl_time(end)]);ylim([40 90]);vline(0,'color','k','linestyle','--','linewidth',2);
text(0.02,43,'CUE onset','HorizontalAlignment','left');title('Distractor RFT (noisy distractors)')
xlabel('Time(s)');ylabel('RFT power (%Baseline)')
text(-0.55,90,'E','Fontsize',20);box off
subplot(4,8,21:24)
h1=fill([dt_time fliplr(dt_time)],[m_RFT_dt_distractor{1}'-se_RFT_dt_distractor{1} fliplr(m_RFT_dt_distractor{1}'+se_RFT_dt_distractor{1})],[0 0 1],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([dt_time fliplr(dt_time)],[m_RFT_dt_distractor{2}'-se_RFT_dt_distractor{2} fliplr(m_RFT_dt_distractor{2}'+se_RFT_dt_distractor{2})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(dt_time,m_RFT_dt_distractor{1},'color',[0 0 1],'linewidth',2,'linestyle','--');
h4=plot(dt_time,m_RFT_dt_distractor{2},'color',[1 0 0],'linewidth',2,'linestyle','--');
xlim([-cl_time(end) 0.25]);ylim([40 90]);vline(0,'color','k','linestyle','--','linewidth',2);
text(-0.02,43,'Discrimination target','HorizontalAlignment','right');legend([h3,h4],{'low target load','high target load'},'location','northwest');box off
xlabel('Time(s)');ax=gca;ax.YTickLabel={};

%SECTION F - Distractor RFT salience effet timecourse
subplot(4,8,25:28)
h1=fill([cl_time fliplr(cl_time)],[m_RFT_distractor{2}'-se_RFT_distractor{2} fliplr(m_RFT_distractor{2}'+se_RFT_distractor{2})],[1 0 0],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([cl_time fliplr(cl_time)],[m_RFT_distractor{4}'-se_RFT_distractor{4} fliplr(m_RFT_distractor{4}'+se_RFT_distractor{4})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(cl_time,m_RFT_distractor{2},'color',[1 0 0],'linewidth',2,'linestyle','--');
h4=plot(cl_time,m_RFT_distractor{4},'color',[1 0 0],'linewidth',2);
xlim([-0.25 cl_time(end)]);ylim([40 90]);vline(0,'color','k','linestyle','--','linewidth',2);
text(0.02,43,'CUE onset','HorizontalAlignment','left');title('Distractor RFT (high load)')
xlabel('Time(s)');ylabel('RFT power (%Baseline)');box off
text(-0.5,90,'F','Fontsize',20)
subplot(4,8,29:32)
h1=fill([dt_time fliplr(dt_time)],[m_RFT_dt_distractor{2}'-se_RFT_dt_distractor{2} fliplr(m_RFT_dt_distractor{2}'+se_RFT_dt_distractor{2})],[1 0 0],'Edgealpha',0);set(h1,'FaceAlpha',0.25);hold on
h2=fill([dt_time fliplr(dt_time)],[m_RFT_dt_distractor{4}'-se_RFT_dt_distractor{4} fliplr(m_RFT_dt_distractor{4}'+se_RFT_dt_distractor{4})],[1 0 0],'Edgealpha',0);set(h2,'FaceAlpha',0.25);
h3=plot(dt_time,m_RFT_dt_distractor{2},'color',[1 0 0],'linewidth',2,'linestyle','--');
h4=plot(dt_time,m_RFT_dt_distractor{4},'color',[1 0 0],'linewidth',2);
xlim([-cl_time(end) 0.25]);ylim([40 90]);vline(0,'color','k','linestyle','--','linewidth',2);
text(-0.02,43,'Discrimination target','HorizontalAlignment','right');legend([h3,h4],{'Noisy distractors','Salient distractors'},'location','northwest')
xlabel('Time(s)');ax=gca;ax.YTickLabel={};box off

%% Calculate correlations with behaviour for paper (mention)
DE_dist_high=(([timecourses.alpha.dist{:,4}]-[timecourses.alpha.dist{:,2}])./([timecourses.alpha.dist{:,4}]+[timecourses.alpha.dist{:,2}]))';
DE_time=timecourses.time+0.35;
t1=dsearchn(DE_time',0.5);
t2=dsearchn(DE_time',1.35);

%Alpha
Alpha_DE_dist_high_tw=nanmean(DE_dist_high(:,t1:t2),2);
[rho,p]=corr(Alpha_DE_dist_high_tw,DE_high,'type','spearman');
disp(['High load Distractor alpha salient vs noisy correlation with behavioral distractor interference with highload: r: ' num2str(rho,2) ' p=' num2str(p,2) ' (spearman)'])

%RFT
RFT_DE_dist_high=(([timecourses.RFT.dist{:,4}]-[timecourses.RFT.dist{:,2}])./([timecourses.RFT.dist{:,4}]+[timecourses.RFT.dist{:,2}]))';
RFT_DE_dist_high_tw=nanmean(RFT_DE_dist_high(:,t1:t2),2);
[rho,p]=corr(RFT_DE_dist_high_tw,DE_high,'type','spearman');
disp(['High load Distractor RFT salient vs noisy correlation with behavioral distractor interference with highload: r: ' num2str(rho,2) ' p=' num2str(p,2) ' (spearman)'])



