function [] = Source_GA()
%
% Generic function to create a source grand average from individual source
% estimates in a folder

%set path
run('/rds/projects/2017/jenseno-01/Tjerk/set_path');


%Alpha
%sel_dir='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/group/Source/Alpha/load_4-3_cl_0.5-1.4/'; %load effect
%sel_dir='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/group/Source/Alpha/load_4-2_cl_0.5-1.1/'; %salience effect

%RFT
%sel_dir='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/group/Source/RFT/RelBL_load_4-3_cl_0.5-1.1/'; %load effect
sel_dir='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/group/Source/RFT/RelBL_load_4-2_cl_0.5-1.1/'; %salience effect

files=dir([sel_dir filesep 'S*.mat']);
cnt=1;
for f=1:length(files)
    sub_nr=str2num(files(f).name(2:3));
    if ~isempty(sub_nr)
        datasets{sub_nr}=[files(f).folder filesep files(f).name];
        cnt=cnt+1;
    end
end

Nsubs=cnt-1;

%load data
for c=1:Nsubs
    disp(['Loading ' datasets{c}])
    tmp=load(datasets{c});
    field=fieldnames(tmp);
    for i=1:length(field)
        src.(field{i}){c}=getfield(tmp,field{i});
    end
end
 
%To equalize marginal differences in the inside field we'll take the first inside field
if isfield(src.(field{1}){1},'inside')
    ins_tmp=src.(field{1}){1}.inside;
else
    ins_tmp=src.(field{1}){1}{1}.inside;
end
for i=1:Nsubs
    for f=1:length(field)
        if isempty(strmatch('plot_data',field{f}))
            if isfield(src.(field{f}){i},'inside')
                src.(field{f}){i}.inside=ins_tmp;
            else
                for f2=1:length(src.(field{f}){i})
                    if isfield(src.(field{f}){i}{f2},'inside')
                        src.(field{f}){i}{f2}.inside=ins_tmp;
                    end
                end
            end
        end
    end
end

if isfield(src.(field{1}){1},'anatomy')
    anatomy=src.(field{1}){1}.anatomy;
else
    anatomy=src.(field{1}){1}{1}.anatomy;
end

plot_data_av=[];

%RFT is still split for 63 and 70 Hz, should be combined before subject averaging
if length(src.(field{1}){1})>1 %two frequencies
    for c=1:Nsubs
        for i=1:length(field)
            src.(field{i}){c}{2}.freq=src.(field{i}){c}{1}.freq; %equalize freqs
            tmp=ft_sourcegrandaverage([],src.(field{i}){c}{1},src.(field{i}){c}{2}); %average freqs
            src.(field{i}){c}=tmp;
        end
    end
end
   
%Average source data
cfg=[];
cfg.keepindividual='yes'; %for subject rejection
for i=1:length(field)
    if isfield(src.(field{i}){1},'pos')
        src_av.(field{i})=ft_sourcegrandaverage(cfg,src.(field{i}){:});
        src_av.(field{i}).anatomy=anatomy;
        src_av.(field{i}).cfg=[];
    else %this may be plottable sensor data
        subfield=fieldnames(src.(field{i}){1});
        for s=1:length(subfield)
            if isfield(src.(field{i}){1}.(subfield{s}),'powspctrm')
                for n=1:Nsubs
                    plot_data.(subfield{s}){n}=getfield(src.(field{i}){n},subfield{s});
                end
                plot_data_av.(subfield{s})=ft_freqgrandaverage(cfg,plot_data.(subfield{s}){:});
                plot_data_av.(subfield{s}).cfg=[];
            end
        end                
    end
end
 
%Save grand average
disp(['Saving ' sel_dir filesep 'Source_GA.mat'])
save([sel_dir filesep 'Source_GA.mat'],'src_av','plot_data_av','-v7.3')
disp('Done')