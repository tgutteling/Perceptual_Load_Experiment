function [] = A2_preICA_cleanup(sub)
%
% To avoid non-convergence of ICA gross artifacts and noisy sensors need to
% be removed.
% This is only a specification of bad trials and sensors
%

%% Set up folders

%proc data location
data_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/';

%Select subject to process if not specified
if nargin<1
    choice_made=0;
    sub_folders=dir([data_folder filesep 'S*']);
    fprintf('Select subject to analyse: \n')
    cnt=1;
    for s=1:size(sub_folders,1)
        if sub_folders(s).isdir
            fprintf(['[' int2str(cnt) '] ' sub_folders(s).name '\n'])
            folder_list(cnt,1)=sub_folders(s);
            cnt=cnt+1;
        end
    end
    while ~choice_made
        choice=input('Enter number: ');
        if choice>size(folder_list,1) || choice<1
            disp('Invalid choice, please try again');
        else
            sub=folder_list(choice).name;
            disp(['Subject ' sub ' selected'])
            choice_made=1;
        end
    end
    choice_made=0;
    fprintf('Select dataset: \n')
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

sub_nr=str2num(sub(2:end));

%load data
disp(['Loading ' data_folder sub filesep sub '_proc.mat'])
load([data_folder sub filesep sub '_proc.mat']);


%get MEG channels
disp('Selecting MEG channels..')
if av_type==1
    MEGchannels=strmatch('MEG',data.meg.label);
    cfg=[];
    cfg.channel=MEGchannels;
    all=ft_selectdata(cfg,data.meg);
else
    MEGchannels=strmatch('MEG',data.meg_dt.label);
    cfg=[];
    cfg.channel=MEGchannels;
    all=ft_selectdata(cfg,data.meg_dt);
end

clear data

%% Visual inspection
disp('Starting visual inspection')
%Sensors and trials will be removed manually, aided by rejectvisual and
%visual inspection per trial. Only gross artifacts (large jumps, excessive
%muscle activity) and excessively noisy sensors will be removed.
ft_rejectvisual([],all);

%note bad chans and trials

%Bad Channels - likely candidates (pre-move): MEG1843 (242) & MEG2122 (261) (these sensors malfunctioned before MEG maintenance)
bad_chans{1}={'MEG1843','MEG2122'};
bad_chans{2}={'MEG1843','MEG2122'};
bad_chans{3}={'MEG1843','MEG2122'};
bad_chans{4}={'MEG1843','MEG2122'};
bad_chans{5}={'MEG1843','MEG2122'};
bad_chans{6}={'MEG1843','MEG2122'};
bad_chans{7}={'MEG1843','MEG2122'};
bad_chans{8}={'MEG1843','MEG2122'};
bad_chans{9}={'MEG1843','MEG2122'};
bad_chans{10}={'MEG1843','MEG2122'};
bad_chans{11}={'MEG1843','MEG2122'};
bad_chans{12}={'MEG1843','MEG2122','MEG0413'};
bad_chans{13}={'MEG1843','MEG2122','MEG0413'};
bad_chans{14}={'MEG1843','MEG2122','MEG0413'};
bad_chans{15}={'MEG1843','MEG2122'};
bad_chans{16}={'MEG1843','MEG2122'};
bad_chans{17}={'MEG0141','MEG1521','MEG1541','MEG1611','MEG1721','MEG2421','MEG2641','MEG1642','MEG1612','MEG1712','MEG1913'};
bad_chans{18}={'MEG1612','MEG0141','MEG1521','MEG1541','MEG1611','MEG2641'};
bad_chans{19}={'MEG2542'};
bad_chans{20}={};
bad_chans{21}={};
bad_chans{22}={'MEG0242'};
bad_chans{23}={};
bad_chans{24}={};
bad_chans{25}={};
bad_chans{26}={};
bad_chans{27}={};
bad_chans{28}={};
bad_chans{29}={};
bad_chans{30}={};
bad_chans{31}={'MEG0413'}; 
bad_chans{32}={'MEG0413'}; 
bad_chans{33}={'MEG0413'}; 
bad_chans{34}={'MEG0413'}; 
bad_chans{35}={'MEG0413'}; 

%% Bad Trials - CUE LOCKED
switch av_type
    case 1 %Cue locked
        bad_trials{1}=[27    59   125   128   141   143   165   209   254   405   407   417   419   422   426   431   433];
        bad_trials{2}=[30    88   113   116   142   144   214   314   454   480   481];
        bad_trials{3}=[110   125   127   129   407   414   421   427   431];
        bad_trials{4}=[140   146   375   393   399   403   412   418   430   439   462   464];
        bad_trials{5}=[6    26    60    73    76    79    82    86   103   121   148   157   176   235   283   344   347   351   365 413   422   483];
        bad_trials{6}=[50    53    67    71    76   331   404   413   417   422   463];
        bad_trials{7}=[169   173   182   185   186   194   195   198   262   335   470   474   479   483   496];
        bad_trials{8}=[36    40    43    72   188   283   327   336   344   346   350   351   373   409   435   460   461   462   463 464   465   466   467   468   490   498   500   502   507   511   512];
        bad_trials{9}=[75    80    87    91   329   383   411   415   422   477   478   479];
        bad_trials{10}=[71   239   250   499   503   507];
        bad_trials{11}=[149   157   166   170   435   439   443   447   451   460];
        bad_trials{12}=[53    79   185   234   235   240   243   251   257   258   364   422   465   491   506];
        bad_trials{13}=[157   302   306   308   310   318   323   325   332];
        bad_trials{14}=[3    89    98   198   233   234   243   244   249   253   302   462];
        bad_trials{15}=[79   156   157   226   227   230   232   235   237   251   409   491];
        bad_trials{16}=[95   108   118   128   226   295   390   409   478   479   483   485   488   489   495   499   500];
        bad_trials{17}=[14	20	45	46	77	114 117 118	129 130	137	157	268	269	271	274	457	458	479	491]; %S17
        bad_trials{18}=[44, 153,154,155,44, 170, 171, 172, 179, 180, 182, 184, 191, 224,225, 280, 281, 302, 303, 304, 335, 344, 349, 350, 351, 362, 363, 384, 387, 393, 403, 404, 407,408, 411, 462, 483,484, 507]; %S18
        bad_trials{19}=[133, 206, 208, 257, 263, 279, 314, 336,  406, 415, 448];%S19
        bad_trials{20}=[53, 66, 69, 100, 105, 131, 183, 209, 279, 283, 291, 292, 309, 313, 335, 358, 362, 368, 387, 388, 391, 405, 406, 407, 413, 429, 435, 465, 468, 469, 470, 471, 491];%S20
        bad_trials{21}=[ 5, 67,  80, 127, 153, 157, 303,  324, 445];%S21
        bad_trials{22}=[19, 23, 25, 26,  119, 120, 156, 183, 184, 209, 273, 274, 283, 309, 341, 342, 465];%S22
        bad_trials{23}=[44, 45, 79, 80, 131, 171, 172, 183, 184, 257, 283, 284, 289, 364, 365, 368, 370, 388, 390, 413, 414, 415, 416, 439, 476, 477, 478, 479, 480, 482, 485, 487, 491, 492];%S23
        bad_trials{24}=[41,  127, 128, 129, 130:140, 142:147,200,205,206, 219, 220, 222, 223,224, 225,226,228:234, 248, 261, 262,263, 276, 326, 335, 336,337, 362,363, 387, 388,395, 400:403, 421,424, 425, 439, 440, 441, 455, 456, 471, 492, 493,504];%S24
        bad_trials{25}=[8, 18, 46, 109, 131, 150, 171, 194, 223, 247, 255, 258:275, 278,300, 309, 310, 351, 359,361:363,369, 373, 385, 386, 387, 394, 417, 456, 458, 459, 476, 479, 480, 501];%sub25
        bad_trials{26}=[151, 152, 160, 213, 477, 483, 506];%S26
        bad_trials{27}=[79, 105, 131, 157, 183, 209, 235, 236, 237, 242, 245, 246, 257, 361, 391, 399, 407, 418, 439, 440,441, 444, 448];
        bad_trials{28}=[79, 183, 235, 413, 439];%S28
        bad_trials{29}=[ 105, 131, 226, 257, 309, 326, 361, 387, 420, 472, 491]; %S29
        bad_trials{30}=[85, 106, 107, 160, 173, 192, 236, 290, 351, 430, 464, 483, 484, 486]; %S30
        bad_trials{31}=[15, 38, 74, 75, 79, 80, 81, 122, 123, 140,141, 183, 293, 306,312, 318, 328,329,331, 335, 358, 365, 372, 387,394, 405, 406, 410, 413, 414, 421, 424, 430,435, 506]; %S31
        bad_trials{32}=[26,126, 131, 221, 252, 257,  283, 361]; %S32
        bad_trials{33}=[52, 55, 80,  106, 157, 158, 165, 171, 175,207,208,209, 210, 215, 218, 226, 233, 239, 246, 247, 262, 277, 312, 325, 335, 336, 343,347, 350, 358, 361, 377, 380, 381, 387, 393, 398, 400,406, 408,409, 411, 413,414, 416, 418, 419, 422, 426, 428, 429, 435, 438, 439, 440, 441, 452, 471, 491, 509, 512]; %S33
        bad_trials{34}=[131, 210,  308, 353, 380, 387, 467, 488, 491]; %S34
        bad_trials{35}=[2, 15,20, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34]; %S35
        
        %% Bad Trials - DT LOCKED
    case 2 %DT Locked        
        bad_trials{1}=[27,113,119,122,124,128,138,141,143,405,414,417,419,421,426,431,433];        
        bad_trials{2}=[116,144,220,401,443,454,469,480,481,493,505,510];        
        bad_trials{3}=[22,42,72,94,110,114,119,127,129,137,139,405,407,413,421,427,433,436,465];
        bad_trials{4}=[132,140,146,148,336,337,393,400,406,412,418];
        bad_trials{5}=[6,10,67,71,76,79,82,86,90,95,103,176,214,283,338,342,343,346,351,355,364,365,413,439];
        bad_trials{6}=[50,57,62,66,67,71,76,173,251,272,361,399,413,417,426,491];
        bad_trials{7}=[166,173,177,182,185,186,189,194,195,198,469,474,478,483,484,488,489,491,496,500];
        bad_trials{8}=[17,22,26,28,32,36,43,72,157,188,283,321,324,325,327,329,331,335,336,337,341,343,346,350,351,373,404,409,411,421,435,460,462,463,466,467,474,477,479,484,486,489,498,507,511,512];
        bad_trials{9}=[75,80,86,95,100,328,329,383,411,415,417,422,429,434,477,478,479];
        bad_trials{10}=[];
        bad_trials{11}=[144,153,157,161,170,439,440,443,444,447,451,455,460];
        bad_trials{12}=[53,185,234,235,238,240,243,258,259,363,364,422,468,505,506];
        bad_trials{13}=[157,302,306,308,310,315,318,322,325,327,332];
        bad_trials{14}=[2,3,89,228,229,233,234,239,243,244,248,249,253,462];
        bad_trials{15}=[79,111,156,183,226,227,230,231,232,235,236,237,238,251,409,443,481,491];
        bad_trials{16}=[95,104,112,118,226,227,254,287,295,305,313,315,324,376,409,435,443,451,478,479,483,488,489,495,500];
        bad_trials{17}=[2,13,14,20,45,46,77,129,136,183,262,268,269,271,335,457,458,479,491];
        bad_trials{18}=[44,155,170,171,172,179,180,182,183,191,192,224,225,279,280,281,302,303,305,335,349,350,351,362,384,387,393,403,407,408,411,462,465,466,483,484];
        bad_trials{19}=[206,208,257,263,314,406,414,447,448];
        bad_trials{20}=[53,66,69,100,183,283,291,292,368,387,388,405,406,407,429,435,438,469,470,471,491];
        bad_trials{21}=5;
        bad_trials{22}=[19,23,25,26,119,120,156,183,184,209,273,274,283,309,341,342,465];
        bad_trials{23}=[79,80,131,171,172,183,283,284,289,364,368,370,388,390,413,414,415,416,439,476,477,478,479,480,482,484,485,487,491];
        bad_trials{24}=[127,128,129,130,131,132,133,134,135,136,137,138,139,140,142,143,144,145,146,148,149,154,177,178,179,200,205,206,219,220,221,222,223,224,225,227,228,229,230,231,232,233,234,248,261,264,275,276,326,335,337,339,342,351,352,355,361,362,363,369,387,395,400,401,403,409,415,421,424,425,439,440,455,456,457,458,465,466,471,492,504,510];
        bad_trials{25}=[8,18,46,62,109,125,131,160,171,194,223,235,247,255,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,309,310,322,340,347,351,352,359,361,362,363,373,385,388,389,393,394,397,412,413,418,422,434,452,456,458,459,476,478,479,480,484,491,495,496,499];
        bad_trials{26}=[150,151,160,213,476,506];
        bad_trials{27}=[53,105,131,157,183,209,236,237,242,245,246,283,308,361,391,397,399,406,407,413,418,439,440,444,448];
        bad_trials{28}=[79,183,235,439];       
        bad_trials{29}=[65,79,105,131,225,226,232,257,309,326,341,361,387,420,422,454,457,472,491];
        bad_trials{30}=[85,106,107,160,172,192,235,236,343,351,404,430,434,447,450,464,470,482,483,484,486,506,510];
        bad_trials{31}=[14,38,51,68,73,75,79,80,122,194,229,232,235,238,293,328,358,365,387,405,410,413,414,421,424,435,457,506];
        bad_trials{32}=[22,26,126,221,257];
        bad_trials{33}=[106,158,165,171,175,182,209,210,233,239,246,247,262,277,312,335,336,350,377,380,381,387,392,398,406,411,413,414,415,416,418,419,422,426,428,438,439,440,441,452,471,491,509,512];
        bad_trials{34}=[131,210,303,308,353,380,387,467,488,491];
        bad_trials{35}=[2,506];
end

%save subject specific.
chan_rej=bad_chans{sub_nr};
trial_rej=bad_trials{sub_nr};

clear bad_trials bad_chans

bad_trials=trial_rej;
bad_chans=chan_rej;

switch av_type
    case 1
        save([data_folder sub filesep sub '_preICA_bad.mat'],'bad_chans','bad_trials');
    case 2
        save([data_folder sub filesep sub '_preICA_bad_dt.mat'],'bad_chans','bad_trials');
end

disp('Done!')
