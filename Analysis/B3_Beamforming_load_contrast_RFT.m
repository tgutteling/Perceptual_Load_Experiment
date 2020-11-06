function [] = B3_Beamforming_load_contrast_RFT(sub_name,bl_window,time_window,contrast)
%
% Here we combine the headmodel with the data to do some source
% localization
%
% Source localizations are estimated using a normalized leadfield and the pre-cue baseline is subtracted, as in the sensor level
% analysis to be able to compare conditions with different noise levels (target high vs low load).
%
%
% Here we want to see localized RFT activity related different load
% comparisons
%
% Time windows should be specified relative to cue onset, e.g. 0.5 - 1.35


disp('Running RFT Beamformer')

%add scripts location
%set path
run('/rds/projects/2017/jenseno-01/Tjerk/set_path');
ft_defaults;

%data folder
proc_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/';

%importing clean data
fprintf(['Loading MEG data: ' sub_name '_all_clean.mat...']);
load([sub_folder sub_name '_all_clean.mat']);disp('Done.')

%get headmodel
load([sub_folder 'headmodel.mat']);

%save grads
grad{1}=data{1}.grad;
grad{2}=data{2}.grad;

%% Use aligned mri to create source model

%getting MNI template sourcemodel
[~,ftdir]=ft_version;
template_dir=[ftdir filesep 'template' filesep 'sourcemodel' filesep];
tmp=load([template_dir 'standard_sourcemodel3d5mm.mat']);
template_grid=tmp.sourcemodel;

%Calculate inverse normalization
cfg                = [];
cfg.grid.warpmni   = 'yes';
cfg.grid.template  = template_grid;
cfg.grid.nonlinear = 'yes'; % use non-linear normalization
cfg.mri            = mri_aligned;
sourcemodel        = ft_prepare_sourcemodel(cfg);

%% Functional data preprocessing

%channel selection
MEG_sens=strmatch('MEG',data{1}.left{1}.meg.label);
grads=[data{1}.left{1}.meg.label(find(str2num(cellfun(@(x) x(end),data{1}.left{1}.meg.label(MEG_sens),'UniformOutput',1))==2)) ;data{1}.left{1}.meg.label(find(str2num(cellfun(@(x) x(end),data{1}.left{1}.meg.label(MEG_sens),'UniformOutput',1))==3))];

%select time window - preCue
disp('Selecting data from baseline time window')
cfg=[];
cfg.toilim = bl_window-0.35;

for config=1:2 %two frequencies
    for l=1:4 %4 load conditions
        bl{config}.left{l}.meg = ft_redefinetrial(cfg,data{config}.left{l}.meg);
        bl{config}.right{l}.meg = ft_redefinetrial(cfg,data{config}.right{l}.meg);
    end
end

%select time window - postCue
disp('Selecting data time window of interest')
cfg=[];
cfg.toilim = time_window-0.35;

for config=1:2 %two frequencies
    for l=1:4 %4 load conditions
        data{config}.left{l}.meg = ft_redefinetrial(cfg,data{config}.left{l}.meg);
        data{config}.right{l}.meg = ft_redefinetrial(cfg,data{config}.right{l}.meg);
    end
end

%create one big dataset for common filter per freq and config
for config=1:2
    dataAll{config}=ft_appenddata([],data{config}.left{contrast(1)}.meg,data{config}.left{contrast(2)}.meg,data{config}.right{contrast(1)}.meg,data{config}.right{contrast(2)}.meg,bl{config}.left{contrast(1)}.meg,bl{config}.left{contrast(2)}.meg,bl{config}.right{contrast(1)}.meg,bl{config}.right{contrast(2)}.meg);
end

RFT_freqs=[63 70];

for freq=1:2
    disp('Calculating CSD');
    cfg=[];
    cfg.channel=grads;
    cfg.method='mtmfft';
    cfg.taper = 'dpss';
    cfg.output = 'powandcsd';
    cfg.keeptrials = 'no';
    cfg.foi=RFT_freqs(freq);
    cfg.tapsmofrq = 2;
    cfg.pad = 2;
    
    
    for config=1:2
        %baseline
        tmp_bl_cond1=ft_appenddata([],bl{config}.left{contrast(1)}.meg,bl{config}.right{contrast(1)}.meg);
        tmp_bl_cond2=ft_appenddata([],bl{config}.left{contrast(2)}.meg,bl{config}.right{contrast(2)}.meg);
        
        cond1.blCSD{config,freq} = ft_freqanalysis(cfg,tmp_bl_cond1);
        cond2.blCSD{config,freq} = ft_freqanalysis(cfg,tmp_bl_cond2);
        
        %post-cue
        cond1.leftCSD{config,freq} = ft_freqanalysis(cfg,data{config}.left{contrast(1)}.meg);
        cond1.rightCSD{config,freq} = ft_freqanalysis(cfg,data{config}.right{contrast(1)}.meg);
        cond2.leftCSD{config,freq} = ft_freqanalysis(cfg,data{config}.left{contrast(2)}.meg);
        cond2.rightCSD{config,freq} = ft_freqanalysis(cfg,data{config}.right{contrast(2)}.meg);
    end
    
    
    %all data to calculate common filter
    for config=1:2
        RFT_all{config,freq} = ft_freqanalysis(cfg,dataAll{config});
    end
end

%% The grad structure (spec. the tra) needs to have the same sensor order as the data
data_labels=data{1}.left{1}.meg.label(MEG_sens);
new_grad = reorder_grad(data_labels,grad{1});

%% Prepare leadfield
disp('Preparing leadfield');
vol_cm = ft_convert_units(vol,'cm');
cfg = [];
cfg.channel=cond1.leftCSD{1,1}.label;
cfg.grad = new_grad;
cfg.headmodel = vol_cm;
cfg.reducerank = 2;
cfg.normalize='yes';

%grid options
cfg.grid.pos=sourcemodel.pos;
cfg.grid.inside=sourcemodel.inside;
cfg.grid.dim=sourcemodel.dim;
cfg.grid.unit = 'cm';

grid = ft_prepare_leadfield(cfg);

%% Calculate common filter
disp('Calculating common filter')
cfg=[];
cfg.method = 'dics';
cfg.grad=grad{1};
cfg.grid=grid;
cfg.headmodel = vol_cm;
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda = '5%';
cfg.dics.keepfilter = 'yes';
cfg.dics.fixedori = 'yes';

for freq=1:2
    cfg.frequency=RFT_freqs(freq);
    for config=1:2
        source_all{config,freq} = ft_sourceanalysis(cfg,RFT_all{config,freq});
    end
end

%% do the sourceanalysis per condition using common filter
disp('Getting to the source')
for freq=1:2
    cfg.frequency=RFT_freqs(freq);
    
    for config=1:2
        
        %set common filter
        cfg.grid.filter=source_all{config,freq}.avg.filter;
        
        %pre-cue
        source_bl.cond1{config,freq} = ft_sourceanalysis(cfg,cond1.blCSD{config,freq});
        source_bl.cond2{config,freq} = ft_sourceanalysis(cfg,cond2.blCSD{config,freq});
        
        %post-cue
        source_left.cond1{config,freq} = ft_sourceanalysis(cfg,cond1.leftCSD{config,freq});
        source_right.cond1{config,freq}= ft_sourceanalysis(cfg,cond1.rightCSD{config,freq});
        source_left.cond2{config,freq} = ft_sourceanalysis(cfg,cond2.leftCSD{config,freq});
        source_right.cond2{config,freq} = ft_sourceanalysis(cfg,cond2.rightCSD{config,freq});
    end
end

%% Source contrasts
disp('Creating contrasts')

%Subtract baseline
for config=1:2
    for freq=1:2
        source_left.cond1{config,freq}.avg.pow=(source_left.cond1{config,freq}.avg.pow-source_bl.cond1{config,freq}.avg.pow);
        source_right.cond1{config,freq}.avg.pow=source_right.cond1{config,freq}.avg.pow-source_bl.cond1{config,freq}.avg.pow;
        source_left.cond2{config,freq}.avg.pow=source_left.cond2{config,freq}.avg.pow-source_bl.cond2{config,freq}.avg.pow;
        source_right.cond2{config,freq}.avg.pow=source_right.cond2{config,freq}.avg.pow-source_bl.cond2{config,freq}.avg.pow;
    end
end

%initizalize contrasts
left_target_contrast{1}=source_left.cond1{1,1};
left_target_contrast{2}=source_left.cond1{1,1};
right_target_contrast{1}=source_right.cond1{1,1};
right_target_contrast{2}=source_right.cond1{1,1};

left_dist_contrast{1}=source_left.cond1{1,1};
left_dist_contrast{2}=source_left.cond1{1,1};
right_dist_contrast{1}=source_right.cond1{1,1};
right_dist_contrast{2}=source_right.cond1{1,1};

%difference
left_target_contrast{1}.avg.pow = (source_left.cond1{1,1}.avg.pow - source_left.cond2{1,1}.avg.pow);
left_target_contrast{2}.avg.pow = (source_left.cond1{2,2}.avg.pow - source_left.cond2{2,2}.avg.pow);
right_target_contrast{1}.avg.pow = (source_right.cond1{2,1}.avg.pow - source_right.cond2{2,1}.avg.pow);
right_target_contrast{2}.avg.pow = (source_right.cond1{1,2}.avg.pow - source_right.cond2{1,2}.avg.pow);

left_dist_contrast{1}.avg.pow = (source_left.cond1{2,1}.avg.pow - source_left.cond2{2,1}.avg.pow);
left_dist_contrast{2}.avg.pow = (source_left.cond1{1,2}.avg.pow - source_left.cond2{1,2}.avg.pow);
right_dist_contrast{1}.avg.pow = (source_right.cond1{1,1}.avg.pow - source_right.cond2{1,1}.avg.pow);
right_dist_contrast{2}.avg.pow = (source_right.cond1{2,2}.avg.pow - source_right.cond2{2,2}.avg.pow);

%average across frequency
att_left_target=ft_sourcegrandaverage([],left_target_contrast{:});
att_right_target=ft_sourcegrandaverage([],right_target_contrast{:});
att_left_dist=ft_sourcegrandaverage([],left_dist_contrast{:});
att_right_dist=ft_sourcegrandaverage([],right_dist_contrast{:});


%% Align to MRI

%restore normalized grid position (e.g. re-normalize)
disp('Restoring normalized grid positions')
att_left_target.pos=template_grid.pos;
att_left_target.dim=template_grid.dim;
att_left_dist.pos=template_grid.pos;
att_left_dist.dim=template_grid.dim;
att_right_target.pos=template_grid.pos;
att_right_target.dim=template_grid.dim;
att_right_dist.pos=template_grid.pos;
att_right_dist.dim=template_grid.dim;

%get SPM template MRI
disp('Retreiving SPM T1 template MRI')
template_mri = ft_read_mri([ftdir filesep 'external' filesep 'spm8' filesep 'templates' filesep 'T1.nii']);

%align source  and T1
disp('Aligning functional and anatomical..')
cfg=[];
cfg.parameter='pow';
att_left_target=ft_sourceinterpolate(cfg, att_left_target,template_mri);
att_right_target=ft_sourceinterpolate(cfg, att_right_target,template_mri);
att_left_dist=ft_sourceinterpolate(cfg, att_left_dist,template_mri);
att_right_dist=ft_sourceinterpolate(cfg, att_right_dist,template_mri);

%% save results
disp('Saving results')
folder_name=['RelBL_load_' int2str(contrast(1)) '-' int2str(contrast(2)) '_cl_' num2str(time_window(1),2) '-' num2str(time_window(2),2)];
if ~exist([proc_folder filesep 'group' filesep 'Source' filesep 'RFT' filesep folder_name])
    mkdir([proc_folder filesep 'group' filesep 'Source' filesep 'RFT' filesep folder_name])
end

save([proc_folder filesep 'group' filesep 'Source' filesep 'RFT' filesep folder_name filesep sub_name '_load_' int2str(contrast(1)) '-' int2str(contrast(2)) '_cl_' num2str(time_window(1),2) '-' num2str(time_window(2),2) '.mat'],'att_left_target','att_right_target','att_left_dist','att_right_dist','-v7.3');
disp('Done!')
end

%AUX function
function [new_grad] = reorder_grad(data_labels,old_grad)

new_grad=old_grad;
new_grad.label=data_labels;
for i=1:length(old_grad.label)
    pos=strmatch(old_grad.label{i},data_labels);
    new_grad.chanori(pos,:)=old_grad.chanori(i,:);
    new_grad.chanpos(pos,:)=old_grad.chanpos(i,:);
    new_grad.chantype(pos,:)=old_grad.chantype(i,:);
    new_grad.chanunit(pos,:)=old_grad.chanunit(i,:);
    new_grad.tra(pos,:)=old_grad.tra(i,:);
end

end

