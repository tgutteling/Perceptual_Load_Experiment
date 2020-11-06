function [] = Source_GA_stats()
%
% Creates final statistical map of source localisation, ready for plotting

%set path
run('/rds/projects/2017/jenseno-01/Tjerk/set_path');

%Set directory containing source GA
%sel_dir='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/group/Source/RFT/RelBL_load_4-3_cl_0.5-1.1/'; %RFT load effect
sel_dir='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/group/Source/RFT/RelBL_load_4-2_cl_0.5-1.1/'; %RFT distractor salience
%sel_dir='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/group/Source/Alpha/load_4-2_cl_0.5-1.1/'; %Alpha distractor salience
%sel_dir='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/group/Source/Alpha/load_4-3_cl_0.5-1.4/'; %Alpha target load

source=load([sel_dir filesep 'Source_GA.mat']);

%% Create clean grand average
bad_subs=[17,18,20,24,29];
field=fieldnames(source.src_av);
for i=1:length(field)
    source.src_av.(field{i}).pow(bad_subs,:)=[];
end

%% Create statistical map
zmap=source.src_av;
for i=1:length(field) 
    sd=nanstd(source.src_av.(field{i}).pow')';
    zmap.(field{i}).pow=mean(source.src_av.(field{i}).pow./sd,1)';   
end

%% Create nifti image         
cfg=[];
cfg.parameter='pow';
for i=1:length(field)
    cfg.filename=[sel_dir field{i} '_z'];
    ft_sourcewrite(cfg,zmap.(field{i}));
end