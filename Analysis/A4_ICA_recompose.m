function [] = A4_ICA_recompose(sub,av_type)
%
% Here the components are loaded and visualy inspected. Components will be
% removed if they are of clear ocular or muscular origin, or exhibit strong 
% non-biological signals.
%

%% Set up folders

%set path
run('/rds/projects/2017/jenseno-01/Tjerk/set_path');
%add scripts location

addpath('/rds/projects/2017/jenseno-01/Tjerk/Load2/scripts/');

%proc data location
data_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/';

if nargin<2
    choice_made=0;
    fprintf('Select averaging type: \n')
    fprintf('[1] Cue-locked\n')
    fprintf('[2] Discrimination target-locked\n')
    while ~choice_made
        choice=input('Enter number: ');
        if choice>2 || choice<1
            disp('Invalid choice, please try again');
        else
            av_type=choice;
            switch choice
                case 1
                    disp('Cue-locked selected')
                case 2
                    disp('Discrimination target-locked selected')
            end
            choice_made=1;
        end
    end
end

%Select subject to process if not specified
if nargin<1
    choice_made=0;
    sub_folders=dir([data_folder filesep 'S*']);
    fprintf('Select subject to analyse: \n')
    cnt=1;
    for s=1:size(sub_folders,1)
        switch av_type
            case 1
                if exist([data_folder sub_folders(s).name filesep sub_folders(s).name '_comp.mat'])>0
                    fprintf(['[' int2str(cnt) '] ' sub_folders(s).name '\n'])
                    sel(cnt)=s;
                    cnt=cnt+1;
                end
            case 2
                if exist([data_folder sub_folders(s).name filesep sub_folders(s).name '_comp_dt.mat'])>0
                    fprintf(['[' int2str(cnt) '] ' sub_folders(s).name '\n'])
                    sel(cnt)=s;
                    cnt=cnt+1;
                end
        end
    end
    while ~choice_made
        choice=input('Enter number: ');
        if choice>(cnt-1) || choice<1
            disp('Invalid choice, please try again');
        else
            sub=sub_folders(sel(choice)).name;
            disp(['Subject ' sub ' selected'])
            choice_made=1;
        end
    end
end

sub_nr=str2num(sub(2:end));

%load data
switch av_type
    case 1
        disp(['Loading ' data_folder sub filesep sub '_comp.mat'])
        load([data_folder sub filesep sub '_comp.mat']);
    case 2
        disp(['Loading ' data_folder sub filesep sub '_comp_dt.mat'])
        load([data_folder sub filesep sub '_comp_dt.mat']);
end

%% create common grad structure    
grad_com=data.grad{1};
grad_com.chanori=(data.grad{1}.chanori+data.grad{2}.chanori)/2;
grad_com.chanpos=(data.grad{1}.chanpos+data.grad{2}.chanpos)/2;
grad_com.coilori=(data.grad{1}.coilori+data.grad{2}.coilori)/2;
grad_com.coilpos=(data.grad{1}.coilpos+data.grad{2}.coilpos)/2;
grad_com.tra=(data.grad{1}.tra+data.grad{2}.tra)/2;
if av_type==1
    data.meg.grad=grad_com;
else
    data.meg_dt.grad=grad_com;
end

%% visual inspection of ICA components
cfg = [];
cfg.channel = 1:5;
cfg.continuous='no';
cfg.viewmode = 'component';
cfg.layout = 'neuromag306mag.lay';
cfg.compscale='local';
ft_databrowser(cfg, data.comp);

%% Note down rejected components

switch av_type
    case 1 %CUE-LOCKED
        comp_rej{1}=[5,6,10,12];              %sub 1
        comp_rej{2}=[1,4,16,17,32];           %sub 2
        comp_rej{3}=[4,10,23,30,34,37];       %sub 3
        comp_rej{4}=[2,3,8,13,20];            %sub 4
        comp_rej{5}=[4,7,9,17,25,31,33,34];   %sub 5
        comp_rej{6}=[4,9,13,15,25,26,30];     %sub 6
        comp_rej{7}=[7,17,18,20];     %       %sub 7
        comp_rej{8}=[1,3,6,14,16,17,20,23,26,28]; %sub 8
        comp_rej{9}=[1,2,3,7,27,28,29,30];    %sub 9
        comp_rej{10}=[1,11,25,28];             %sub 10
        comp_rej{11}=[1,22,23,24,30];          %sub 11
        comp_rej{12}=[3,4,11,12,17,20,21,22,24]; %sub 12
        comp_rej{13}=[25,29,36,37,38,40];       %sub 13
        comp_rej{14}=[1,2,5,7,8,9,10,12,14,15,34]; %sub 14
        comp_rej{15}=[1,2,4,5,6,8,12,15,17,19,31,33,34,50]; %sub 15
        comp_rej{16}=[8,10,13,16,19,23,25,28,36,37,39,43]; %sub 16
        comp_rej{17}=[1,4,15,16,17,19,22]; %sub 17
        comp_rej{18}=[3,8,13,18,22,24,29,32,35,37,46,48]; %sub 18
        comp_rej{19}=[6,16,22,24,31]; %sub19
        comp_rej{20}=[2,12,13,15,18,19,20,21]; %sub20
        comp_rej{21}=[19,26,27,28,33,34,35,36,37,38,43]; %sub21
        comp_rej{22}=[9,28,29,33,34,35]; %sub22
        comp_rej{23}=[1,6,8,9,10,12,13,22]; %sub23
        comp_rej{24}=[1,2,3,4,5,6,10,13,15,17,19,20,21,23,24,25,28,36]; %sub24
        comp_rej{25}=[1,2,6,10,15,16,17,18,19,21,26,32,35,37,38,40]; %sub 25                            
        comp_rej{26}=[13,15,19,23,28,31,32,35];%sub 26
        comp_rej{27}=[2,5,6,19,20,21,23,24,38];%sub 27
        comp_rej{28}=[1,2,18,29,41];%sub 28
        comp_rej{29}=[13,21,23,24,27,28,45];%sub 29
        comp_rej{30}=[1:6,9,15,22,25,26];%sub 30
        comp_rej{31}=[1:5,8:10,21,31,41,45]; %sub 31
        comp_rej{32}=[9,24,25,29,32]; %sub32
        comp_rej{33}=[11,12,13,16,17,20,21,22]; %sub33
        comp_rej{34}=[1,3,6,7,11]; %sub34
        comp_rej{35}=[1,2,3,6,8,9,13,15,17]; %sub35
        
        %actually remove components
        cfg=[];
        cfg.component=comp_rej{sub_nr};
        data.meg=ft_rejectcomponent(cfg,data.comp,data.meg); 
        
        %save updated common grad structure
        data.grad_common=data.meg.grad;

        %add back aux channels
        data.meg=ft_appenddata([],data.meg,data.aux);

        %since merged, we can now remove the aux again
        data=rmfield(data,'aux');

        %since we recomposed, remove the ICA components
        data=rmfield(data,'comp');
        
    case 2 %DT LOCKED
        comp_rej{1}=[5,7,10,12];
        comp_rej{2}=[1,4,11,20,21,23,27,32];
        comp_rej{3}=[4,10,21,22,30,34,37,38,41];
        comp_rej{4}=[2,3,9,16,17,23];
        comp_rej{5}=[1,4,5,10,24,35,40];
        comp_rej{6}=[3,9,14,17,18,25,27];    
        comp_rej{7}=[8,15,17,18,27];
        comp_rej{8}=[1,4,6,14,15,16,17,20,33,23];
        comp_rej{9}=[1,2,4,5,7,23,29];
        comp_rej{10}=[1,9,11,16,24,25];
        comp_rej{11}=[8,21,23,24,29];
        comp_rej{12}=[4,5,11,12,17,21,22,24];
        comp_rej{13}=[26,31,33,38,36,40];
        comp_rej{14}=[1,2,8,9,10,11,14];
        comp_rej{15}=[1,2,3,5,6,8,10,15,21,30,32];
        comp_rej{16}=[6,8,14,18,20,23,25,26,33,35,36,37];
        comp_rej{17}=[1,5,13,14,15,20,21];
        comp_rej{18}=[4,10,14,18,20,33,47];
        comp_rej{19}=[9,16,19,22,25,30];
        comp_rej{20}=[5,9,11,14,17,18,20,21];
        comp_rej{21}=[21,22,23,27,35,36,37,40,42];
        comp_rej{22}=[8,27,29,33,35];
        comp_rej{23}=[1,5,9,12,14,26];
        comp_rej{24}=[1:5,7,8,13,14,15,17:21];
        comp_rej{25}=[1,2,6,7,14:17,19,25,28,30,33,34];
        comp_rej{26}=[10,12,13,24,28,29,31,35,36];
        comp_rej{27}=[2,5,7,15,16,22,26,28,37];
        comp_rej{28}=[1,3,16,29];
        comp_rej{29}=[15,19,20,22,26,27];
        comp_rej{30}=[1:7,16,18,26];
        comp_rej{31}=[1:5,7:10,19,20,30,38,43];
        comp_rej{32}=[10,18,26,29,30,32];
        comp_rej{33}=[13,16,17,18,21,22,23];
        comp_rej{34}=[1,3,5,7,10,37];
        comp_rej{35}=[1:9,17,19];
        
        cfg=[];
        cfg.component=comp_rej{sub_nr};
        data.meg_dt=ft_rejectcomponent(cfg,data.comp,data.meg_dt);
        
        %save updated common grad structure
        data.grad_common=data.meg_dt.grad;
        
        %add back aux channels
        data.meg_dt=ft_appenddata([],data.meg_dt,data.aux);

        %since merged, we can now remove the aux again
        data=rmfield(data,'aux');

        %since we recomposed, remove the ICA components
        data=rmfield(data,'comp');        
end

%% Save
fprintf(['Saving data for ' sub '...'])
switch av_type
    case 1
        save([data_folder sub filesep sub '_icclean.mat'],'data','-v7.3');
    case 2
        save([data_folder sub filesep sub '_icclean_dt.mat'],'data','-v7.3');
end
fprintf('..Done\n')



