function [] = B4_Beamforming_load_contrast_alpha(sub_name,time_window,contrast,type)
%
% Here we combine the headmodel with the data to do some source
% localization
%
% Here we want to see localized Alpha activity related different load
% comparisons
%
% Time windows should be specified relative to cue onset, e.g. 0.5 - 1.35

disp('Running Alpha Beamformer (load-contrasts)')

%add scripts location
%set path
run('/rds/projects/2017/jenseno-01/Tjerk/set_path');
ft_defaults;

%data folder
proc_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/';

%select subject
if nargin<1
    sub_folders=dir([proc_folder filesep 'S*']);
    cnt=1;
    for s=1:size(sub_folders,1)
        if exist([proc_folder sub_folders(s).name filesep 'headmodel.mat'])>0
            datasets{cnt}=[proc_folder sub_folders(s).name filesep];
            subject_name{cnt}=sub_folders(s).name;
            cnt=cnt+1;
        end
    end
    choice_made=0;
    
    fprintf('Select Subject to analyse: \n')
    for s=1:cnt-1
        fprintf(['[' int2str(s) '] ' subject_name{s} '\n'])
    end
    while ~choice_made
        choice=input('Enter number: ');
        if choice>size(datasets,2) || choice<1
            disp('Invalid choice, please try again');
        else
            sub_folder=datasets{choice};
            sub_name=subject_name{choice};
            disp(['Subject ' sub_name ' selected'])
            choice_made=1;
        end
    end
else
    sub_folder=[proc_folder sub_name filesep];
end

%importing clean data
switch type
    case 1 %CL
        fprintf(['Loading MEG data: ' sub_name '_all_clean.mat...']);
        load([sub_folder sub_name '_all_clean.mat']);disp('Done.')
    case 2 %DT
        fprintf(['Loading MEG data: ' sub_name '_all_clean_dt.mat...']);
        load([sub_folder sub_name '_all_clean_dt.mat']);disp('Done.')
end

%get headmodel
load([sub_folder 'headmodel.mat']);

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
if type==1
    MEG_sens=strmatch('MEG',data{1}.left{1}.meg.label);    
    grads=[data{1}.left{1}.meg.label(find(str2num(cellfun(@(x) x(end),data{1}.left{1}.meg.label(MEG_sens),'UniformOutput',1))==2)) ;data{1}.left{1}.meg.label(find(str2num(cellfun(@(x) x(end),data{1}.left{1}.meg.label(MEG_sens),'UniformOutput',1))==3))];
else
    MEG_sens=strmatch('MEG',data{1}.left{1}.meg_dt.label);    
    grads=[data{1}.left{1}.meg_dt.label(find(str2num(cellfun(@(x) x(end),data{1}.left{1}.meg_dt.label(MEG_sens),'UniformOutput',1))==2)) ;data{1}.left{1}.meg_dt.label(find(str2num(cellfun(@(x) x(end),data{1}.left{1}.meg_dt.label(MEG_sens),'UniformOutput',1))==3))];
end

%select time window
disp('Selecting data time window of interest')
cfg=[];
cfg.toilim = time_window-0.35;
for config=1:2 %two frequencies
    for l=1:4 %4 load conditions
        if type==1
            data{config}.left{l}.meg = ft_redefinetrial(cfg,data{config}.left{l}.meg);
            data{config}.right{l}.meg = ft_redefinetrial(cfg,data{config}.right{l}.meg);
        else
            data{config}.left{l}.meg_dt = ft_redefinetrial(cfg,data{config}.left{l}.meg_dt);
            data{config}.right{l}.meg_dt = ft_redefinetrial(cfg,data{config}.right{l}.meg_dt);
        end
    end
end

%collapse over RFT config
if type==1
    cond1.left = ft_appenddata([],data{1}.left{contrast(1)}.meg,data{2}.left{contrast(1)}.meg);
    cond1.right = ft_appenddata([],data{1}.right{contrast(1)}.meg,data{2}.right{contrast(1)}.meg);
    cond2.left = ft_appenddata([],data{1}.left{contrast(2)}.meg,data{2}.left{contrast(2)}.meg);
    cond2.right = ft_appenddata([],data{1}.right{contrast(2)}.meg,data{2}.right{contrast(2)}.meg);
else
    cond1.left = ft_appenddata([],data{1}.left{contrast(1)}.meg_dt,data{2}.left{contrast(1)}.meg_dt);
    cond1.right = ft_appenddata([],data{1}.right{contrast(1)}.meg_dt,data{2}.right{contrast(1)}.meg_dt);
    cond2.left = ft_appenddata([],data{1}.left{contrast(2)}.meg_dt,data{2}.left{contrast(2)}.meg_dt);
    cond2.right = ft_appenddata([],data{1}.right{contrast(2)}.meg_dt,data{2}.right{contrast(2)}.meg_dt);
end

%create one big dataset
dataAll=ft_appenddata([],cond1.left,cond1.right,cond2.left,cond2.right);

%correct foi to fit integer nr of cycles
target_foi=10.5;
corr_foi=round((time_window(2)-time_window(1))/(1/target_foi))/(time_window(2)-time_window(1));

disp('Calculating CSD');
cfg=[];
cfg.channel=grads;
cfg.method='mtmfft';
cfg.taper = 'dpss';
cfg.output = 'powandcsd';
cfg.keeptrials = 'no';
cfg.foi=corr_foi;
cfg.tapsmofrq = 3;
cfg.pad = 'nextpow2';

cond1.leftCSD = ft_freqanalysis(cfg,cond1.left);
cond1.rightCSD = ft_freqanalysis(cfg,cond1.right);
cond2.leftCSD = ft_freqanalysis(cfg,cond2.left);
cond2.rightCSD = ft_freqanalysis(cfg,cond2.right);

%all data to calculate common filter
Alpha_all = ft_freqanalysis(cfg,dataAll);

%collapse over frequency
cfg=[];
cfg.avgoverfreq='yes';
cond1.leftCSD=ft_selectdata(cfg,cond1.leftCSD);
cond1.rightCSD=ft_selectdata(cfg,cond1.rightCSD);
cond2.leftCSD=ft_selectdata(cfg,cond2.leftCSD);
cond2.rightCSD=ft_selectdata(cfg,cond2.rightCSD);
Alpha_all=ft_selectdata(cfg,Alpha_all);

%% The grad structure (spec. the tra) needs to have the same sensor order as the data
if type==1
    data_labels=data{1}.left{1}.meg.label(MEG_sens);
else
    data_labels=data{1}.left{1}.meg_dt.label(MEG_sens);
end
new_grad = reorder_grad(data_labels,grad_common);

%% Prepare leadfield
disp('Preparing leadfield');
vol_cm = ft_convert_units(vol,'cm');
cfg = [];
cfg.channel=cond1.leftCSD.label;
cfg.grad = new_grad;
cfg.headmodel = vol_cm;
cfg.reducerank = 2;
cfg.normalize='no';

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
cfg.grad=new_grad;
cfg.grid=grid;
cfg.headmodel = vol_cm;
cfg.dics.lambda = '5%';
cfg.dics.keepfilter = 'yes';
cfg.dics.fixedori = 'yes';
cfg.frequency=cond1.leftCSD.freq;

source_all = ft_sourceanalysis(cfg,Alpha_all);

%% do the sourceanalysis per condition using common filter
disp('Getting to the source')
cfg.grid.filter=source_all.avg.filter;    

source_left.cond1 = ft_sourceanalysis(cfg,cond1.leftCSD);
source_right.cond1 = ft_sourceanalysis(cfg,cond1.rightCSD);
source_left.cond2 = ft_sourceanalysis(cfg,cond2.leftCSD);
source_right.cond2 = ft_sourceanalysis(cfg,cond2.rightCSD);

%% Source contrasts
disp('Creating contrasts')

%initizalize contrasts
left_contrast=source_left.cond1;
right_contrast=source_right.cond1;

left_contrast_index=source_left.cond1;
right_contrast_index=source_right.cond1;

%simple difference
left_contrast.avg.pow = source_left.cond1.avg.pow - source_left.cond2.avg.pow;
right_contrast.avg.pow = source_right.cond1.avg.pow - source_left.cond2.avg.pow;

%difference over sum
left_contrast_index.avg.pow = (source_left.cond1.avg.pow - source_left.cond2.avg.pow) ./ ((source_left.cond1.avg.pow+source_left.cond2.avg.pow)*0.5);
right_contrast_index.avg.pow = (source_right.cond1.avg.pow - source_left.cond2.avg.pow) ./ ((source_left.cond1.avg.pow+source_right.cond2.avg.pow)*0.5);

%% Align to MRI

%restore normalized grid position (e.g. re-normalize)
disp('Restoring normalized grid positions')
left_contrast.pos=template_grid.pos;
left_contrast.dim=template_grid.dim;
right_contrast.pos=template_grid.pos;
right_contrast.dim=template_grid.dim;

left_contrast_index.pos=template_grid.pos;
left_contrast_index.dim=template_grid.dim;
right_contrast_index.pos=template_grid.pos;
right_contrast_index.dim=template_grid.dim;

%get SPM template MRI
disp('Retreiving SPM T1 template MRI')
template_mri = ft_read_mri([ftdir filesep 'external' filesep 'spm8' filesep 'templates' filesep 'T1.nii']);

%align source  and T1
disp('Aligning functional and anatomical..')
cfg=[];
cfg.parameter='pow';
left_contrast=ft_sourceinterpolate(cfg, left_contrast,template_mri);
right_contrast=ft_sourceinterpolate(cfg, right_contrast,template_mri);
left_contrast_index=ft_sourceinterpolate(cfg, left_contrast_index,template_mri);
right_contrast_index=ft_sourceinterpolate(cfg, right_contrast_index,template_mri);

%% save results
disp('Saving results')
type_names={'cl','dt'};
folder_name=['load_' int2str(contrast(1)) '-' int2str(contrast(2)) '_' type_names{type} '_' num2str(time_window(1),2) '-' num2str(time_window(2),2)];
if ~exist([proc_folder filesep 'group' filesep 'Source' filesep 'Alpha' filesep folder_name])
    mkdir([proc_folder filesep 'group' filesep 'Source' filesep 'Alpha' filesep folder_name])
end

save([proc_folder filesep 'group' filesep 'Source' filesep 'Alpha' filesep folder_name filesep sub_name '_load_' int2str(contrast(1)) '-' int2str(contrast(2)) '_' type_names{type} '_' num2str(time_window(1),2) '-' num2str(time_window(2),2) '.mat'],'left_contrast','right_contrast','left_contrast_index','right_contrast_index','-v7.3');
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
