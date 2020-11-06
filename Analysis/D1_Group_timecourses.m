function [] = D1_Group_timecourses()
%
% Final processing on the timecourse alpha and RFT data. 
% Here 

%use only correct trials?
correct_only=1;
%set_path
run('/rds/projects/2017/jenseno-01/Tjerk/set_path');

%folders
group_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/group/';


%% Load data
disp('Loading data')
%RFT
if correct_only
    tmp=load([group_folder 'ROI_RFT_correct_only.mat']);
else
    tmp=load([group_folder 'ROI_RFT.mat']);
end

RFT=tmp.RFT;
RFT_time=tmp.time;

%RFT_dt
if correct_only
    tmp=load([group_folder 'ROI_RFT_dt_correct_only.mat']);
else
    tmp=load([group_folder 'ROI_RFT_dt.mat']);
end

RFT_dt=tmp.RFT_dt;
RFT_dt_time=tmp.time;

%alpha
if correct_only
    tmp=load([group_folder 'ROI_alpha_correct_only.mat']);
else
    tmp=load([group_folder 'ROI_alpha.mat']);
end

alpha=tmp.alpha;
alpha_time=tmp.time;

%alpha_dt
if correct_only
    tmp=load([group_folder 'ROI_alpha_dt_correct_only.mat']);
else
    tmp=load([group_folder 'ROI_alpha_dt.mat']);
end

alpha_dt=tmp.alpha;
alpha_dt_time=tmp.time;

%% Results per load condition

%subject rejection 
bad_subs=[17,18,20,24,29]; 

d=1;
for n=1:length(RFT)
    
    if ~ismember(n,bad_subs)
        disp(['Processing subject ' int2str(n)])
        cur_RFT=RFT{n}.ind;
        cur_RFT_dt=RFT_dt{n}.ind;
        cur_alpha=alpha{n};
        cur_alpha_dt=alpha_dt{n};
        
        %%%%%%%
        % AMI %
        %%%%%%%%
        for l=1:4
            ami_alpha{l}(d,:)=((cur_alpha.target_left{l}-cur_alpha.dist_left{l})./(cur_alpha.target_left{l}+cur_alpha.dist_left{l})+(cur_alpha.target_right{l}-cur_alpha.dist_right{l})./(cur_alpha.target_right{l}+cur_alpha.dist_right{l}))/2;
            ami_alpha_dt{l}(d,:)=((cur_alpha_dt.target_left{l}-cur_alpha_dt.dist_left{l})./(cur_alpha_dt.target_left{l}+cur_alpha_dt.dist_left{l})+(cur_alpha_dt.target_right{l}-cur_alpha_dt.dist_right{l})./(cur_alpha_dt.target_right{l}+cur_alpha_dt.dist_right{l}))/2;
        end
        
        %%%%%%%
        % LMI %
        %%%%%%%
        
        %RFT
        for freq=1:2
            RFT_target_LMI_noisy_dist(d,freq,:)=mean([(cur_RFT.target_left{2,freq}-cur_RFT.target_left{1,freq})./(cur_RFT.target_left{2,freq}+cur_RFT.target_left{1,freq}) (cur_RFT.target_right{2,freq}-cur_RFT.target_right{1,freq})./(cur_RFT.target_right{2,freq}+cur_RFT.target_right{1,freq})],2);
            RFT_target_LMI_salient_dist(d,freq,:)=mean([(cur_RFT.target_left{4,freq}-cur_RFT.target_left{3,freq})./(cur_RFT.target_left{4,freq}+cur_RFT.target_left{3,freq}) (cur_RFT.target_right{4,freq}-cur_RFT.target_right{3,freq})./(cur_RFT.target_right{4,freq}+cur_RFT.target_right{3,freq})],2);
            RFT_dist_LMI_noisy_dist(d,freq,:)=mean([(cur_RFT.dist_left{2,freq}-cur_RFT.dist_left{1,freq})./(cur_RFT.dist_left{2,freq}+cur_RFT.dist_left{1,freq}) (cur_RFT.dist_right{2,freq}-cur_RFT.dist_right{1,freq})./(cur_RFT.dist_right{2,freq}+cur_RFT.dist_right{1,freq})],2);
            RFT_dist_LMI_salient_dist(d,freq,:)=mean([(cur_RFT.dist_left{4,freq}-cur_RFT.dist_left{3,freq})./(cur_RFT.dist_left{4,freq}+cur_RFT.dist_left{3,freq}) (cur_RFT.dist_right{4,freq}-cur_RFT.dist_right{3,freq})./(cur_RFT.dist_right{4,freq}+cur_RFT.dist_right{3,freq})],2);
            
            RFT_dt_target_LMI_noisy_dist(d,freq,:)=mean([(cur_RFT_dt.target_left{2,freq}-cur_RFT_dt.target_left{1,freq})./(cur_RFT_dt.target_left{2,freq}+cur_RFT_dt.target_left{1,freq}) (cur_RFT_dt.target_right{2,freq}-cur_RFT_dt.target_right{1,freq})./(cur_RFT_dt.target_right{2,freq}+cur_RFT_dt.target_right{1,freq})],2);
            RFT_dt_target_LMI_salient_dist(d,freq,:)=mean([(cur_RFT_dt.target_left{4,freq}-cur_RFT_dt.target_left{3,freq})./(cur_RFT_dt.target_left{4,freq}+cur_RFT_dt.target_left{3,freq}) (cur_RFT_dt.target_right{4,freq}-cur_RFT_dt.target_right{3,freq})./(cur_RFT_dt.target_right{4,freq}+cur_RFT_dt.target_right{3,freq})],2);
            RFT_dt_dist_LMI_noisy_dist(d,freq,:)=mean([(cur_RFT_dt.dist_left{2,freq}-cur_RFT_dt.dist_left{1,freq})./(cur_RFT_dt.dist_left{2,freq}+cur_RFT_dt.dist_left{1,freq}) (cur_RFT_dt.dist_right{2,freq}-cur_RFT_dt.dist_right{1,freq})./(cur_RFT_dt.dist_right{2,freq}+cur_RFT_dt.dist_right{1,freq})],2);
            RFT_dt_dist_LMI_salient_dist(d,freq,:)=mean([(cur_RFT_dt.dist_left{4,freq}-cur_RFT_dt.dist_left{3,freq})./(cur_RFT_dt.dist_left{4,freq}+cur_RFT_dt.dist_left{3,freq}) (cur_RFT_dt.dist_right{4,freq}-cur_RFT_dt.dist_right{3,freq})./(cur_RFT_dt.dist_right{4,freq}+cur_RFT_dt.dist_right{3,freq})],2);
        end
        
        %alpha
        Alpha_target_LMI_noisy_dist(d,:)=mean([(cur_alpha.target_left{2}-cur_alpha.target_left{1})./(cur_alpha.target_left{2}+cur_alpha.target_left{1}) (cur_alpha.target_right{2}-cur_alpha.target_right{1})./(cur_alpha.target_right{2}+cur_alpha.target_right{1})],2);
        Alpha_target_LMI_salient_dist(d,:)=mean([(cur_alpha.target_left{4}-cur_alpha.target_left{3})./(cur_alpha.target_left{4}+cur_alpha.target_left{3}) (cur_alpha.target_right{4}-cur_alpha.target_right{3})./(cur_alpha.target_right{4}+cur_alpha.target_right{3})],2);
        Alpha_dist_LMI_noisy_dist(d,:)=mean([(cur_alpha.dist_left{2}-cur_alpha.dist_left{1})./(cur_alpha.dist_left{2}+cur_alpha.dist_left{1}) (cur_alpha.dist_right{2}-cur_alpha.dist_right{1})./(cur_alpha.dist_right{2}+cur_alpha.dist_right{1})],2);
        Alpha_dist_LMI_salient_dist(d,:)=mean([(cur_alpha.dist_left{4}-cur_alpha.dist_left{3})./(cur_alpha.dist_left{4}+cur_alpha.dist_left{3}) (cur_alpha.dist_right{4}-cur_alpha.dist_right{3})./(cur_alpha.dist_right{4}+cur_alpha.dist_right{3})],2);
        
        %alpha_dt
        Alpha_dt_target_LMI_noisy_dist(d,:)=mean([(cur_alpha_dt.target_left{2}-cur_alpha_dt.target_left{1})./(cur_alpha_dt.target_left{2}+cur_alpha_dt.target_left{1}) (cur_alpha_dt.target_right{2}-cur_alpha_dt.target_right{1})./(cur_alpha_dt.target_right{2}+cur_alpha_dt.target_right{1})],2);
        Alpha_dt_target_LMI_salient_dist(d,:)=mean([(cur_alpha_dt.target_left{4}-cur_alpha_dt.target_left{3})./(cur_alpha_dt.target_left{4}+cur_alpha_dt.target_left{3}) (cur_alpha_dt.target_right{4}-cur_alpha_dt.target_right{3})./(cur_alpha_dt.target_right{4}+cur_alpha_dt.target_right{3})],2);
        Alpha_dt_dist_LMI_noisy_dist(d,:)=mean([(cur_alpha_dt.dist_left{2}-cur_alpha_dt.dist_left{1})./(cur_alpha_dt.dist_left{2}+cur_alpha_dt.dist_left{1}) (cur_alpha_dt.dist_right{2}-cur_alpha_dt.dist_right{1})./(cur_alpha_dt.dist_right{2}+cur_alpha_dt.dist_right{1})],2);
        Alpha_dt_dist_LMI_salient_dist(d,:)=mean([(cur_alpha_dt.dist_left{4}-cur_alpha_dt.dist_left{3})./(cur_alpha_dt.dist_left{4}+cur_alpha_dt.dist_left{3}) (cur_alpha_dt.dist_right{4}-cur_alpha_dt.dist_right{3})./(cur_alpha_dt.dist_right{4}+cur_alpha_dt.dist_right{3})],2);
        
        %%%%%%%
        % DE %
        %%%%%%%
        
        %Distractor effect (salient distractor - noisy distractor) / sum
        
        %RFT
        for freq=1:2
            RFT_target_DE_low_load(d,freq,:)=mean([(cur_RFT.target_left{3,freq}-cur_RFT.target_left{1,freq})./(cur_RFT.target_left{3,freq}+cur_RFT.target_left{1,freq}) (cur_RFT.target_right{3,freq}-cur_RFT.target_right{1,freq})./(cur_RFT.target_right{3,freq}+cur_RFT.target_right{1,freq})],2);
            RFT_target_DE_high_load(d,freq,:)=mean([(cur_RFT.target_left{4,freq}-cur_RFT.target_left{2,freq})./(cur_RFT.target_left{4,freq}+cur_RFT.target_left{2,freq}) (cur_RFT.target_right{4,freq}-cur_RFT.target_right{2,freq})./(cur_RFT.target_right{4,freq}+cur_RFT.target_right{2,freq})],2);
            RFT_dt_target_DE_low_load(d,freq,:)=mean([(cur_RFT_dt.target_left{3,freq}-cur_RFT_dt.target_left{1,freq})./(cur_RFT_dt.target_left{3,freq}+cur_RFT_dt.target_left{1,freq}) (cur_RFT_dt.target_right{3,freq}-cur_RFT_dt.target_right{1,freq})./(cur_RFT_dt.target_right{3,freq}+cur_RFT_dt.target_right{1,freq})],2);
            RFT_dt_target_DE_high_load(d,freq,:)=mean([(cur_RFT_dt.target_left{4,freq}-cur_RFT_dt.target_left{2,freq})./(cur_RFT_dt.target_left{4,freq}+cur_RFT_dt.target_left{2,freq}) (cur_RFT_dt.target_right{4,freq}-cur_RFT_dt.target_right{2,freq})./(cur_RFT_dt.target_right{4,freq}+cur_RFT_dt.target_right{2,freq})],2);
        end
        
        %alpha
        Alpha_target_DE_low_load(d,:)=mean([(cur_alpha.target_left{3}-cur_alpha.target_left{1})./(cur_alpha.target_left{3}+cur_alpha.target_left{1}) (cur_alpha.target_right{3}-cur_alpha.target_right{1})./(cur_alpha.target_right{3}+cur_alpha.target_right{1})],2);
        Alpha_target_DE_high_load(d,:)=mean([(cur_alpha.target_left{4}-cur_alpha.target_left{2})./(cur_alpha.target_left{4}+cur_alpha.target_left{2}) (cur_alpha.target_right{4}-cur_alpha.target_right{2})./(cur_alpha.target_right{4}+cur_alpha.target_right{2})],2);
        Alpha_dt_target_DE_low_load(d,:)=mean([(cur_alpha_dt.target_left{3}-cur_alpha_dt.target_left{1})./(cur_alpha_dt.target_left{3}+cur_alpha_dt.target_left{1}) (cur_alpha_dt.target_right{3}-cur_alpha_dt.target_right{1})./(cur_alpha_dt.target_right{3}+cur_alpha_dt.target_right{1})],2);
        Alpha_dt_target_DE_high_load(d,:)=mean([(cur_alpha_dt.target_left{4}-cur_alpha_dt.target_left{2})./(cur_alpha_dt.target_left{4}+cur_alpha_dt.target_left{2}) (cur_alpha_dt.target_right{4}-cur_alpha_dt.target_right{2})./(cur_alpha_dt.target_right{4}+cur_alpha_dt.target_right{2})],2);
                
        Alpha_dist_DE_low_load(d,:)=mean([(cur_alpha.dist_left{3}-cur_alpha.dist_left{1})./(cur_alpha.dist_left{3}+cur_alpha.dist_left{1}) (cur_alpha.dist_right{3}-cur_alpha.dist_right{1})./(cur_alpha.dist_right{3}+cur_alpha.dist_right{1})],2);
        Alpha_dt_dist_DE_low_load(d,:)=mean([(cur_alpha_dt.dist_left{3}-cur_alpha_dt.dist_left{1})./(cur_alpha_dt.dist_left{3}+cur_alpha_dt.dist_left{1}) (cur_alpha_dt.dist_right{3}-cur_alpha_dt.dist_right{1})./(cur_alpha_dt.dist_right{3}+cur_alpha_dt.dist_right{1})],2);
        Alpha_dist_DE_high_load(d,:)=mean([(cur_alpha.dist_left{4}-cur_alpha.dist_left{2})./(cur_alpha.dist_left{4}+cur_alpha.dist_left{2}) (cur_alpha.dist_right{4}-cur_alpha.dist_right{2})./(cur_alpha.dist_right{4}+cur_alpha.dist_right{2})],2);
        Alpha_dt_dist_DE_high_load(d,:)=mean([(cur_alpha_dt.dist_left{4}-cur_alpha_dt.dist_left{2})./(cur_alpha_dt.dist_left{4}+cur_alpha_dt.dist_left{2}) (cur_alpha_dt.dist_right{4}-cur_alpha_dt.dist_right{2})./(cur_alpha_dt.dist_right{4}+cur_alpha_dt.dist_right{2})],2);
        
        %%%%%%%%%%%%%%%
        % Timecourses %
        %%%%%%%%%%%%%%%
        
        %RFT
        %Get pre-RFT baseline
        t1=find(~isnan(cur_RFT.target_left{1,1}),1,'first');%first datapoint
        t2=dsearchn(RFT_time',-1.4); %before RFT onset
        all_left=cell(1,2);
        all_right=cell(1,2);
        for config=1:2
            for l=1:4
                all_left{config}=[all_left{config} ; cur_RFT.target_left{l,config}(t1:t2) ; cur_RFT.dist_left{l,config}(t1:t2)];
                all_right{config}=[all_right{config} ; cur_RFT.target_right{l,config}(t1:t2) ; cur_RFT.dist_right{l,config}(t1:t2)];
            end
        end
        baseline_left=mean([all_left{1} ; all_left{2}]);
        baseline_right=mean([all_right{1} ; all_right{2}]);
        
        %remove baseline
        for freq=1:2
            for l=1:4
                cur_RFT.target_left{l,freq} = cur_RFT.target_left{l,freq}-baseline_left;
                cur_RFT.target_right{l,freq} = cur_RFT.target_right{l,freq}-baseline_right;
                cur_RFT.dist_left{l,freq} = cur_RFT.dist_left{l,freq}-baseline_left;
                cur_RFT.dist_right{l,freq} = cur_RFT.dist_right{l,freq}-baseline_right;
                
                cur_RFT_dt.target_left{l,freq} = cur_RFT_dt.target_left{l,freq}-baseline_left;
                cur_RFT_dt.target_right{l,freq} = cur_RFT_dt.target_right{l,freq}-baseline_right;
                cur_RFT_dt.dist_left{l,freq} = cur_RFT_dt.dist_left{l,freq}-baseline_left;
                cur_RFT_dt.dist_right{l,freq} = cur_RFT_dt.dist_right{l,freq}-baseline_right;
            end
        end
        
        %Determine scaling, use pre-cue interval
        t1=dsearchn(RFT_time',-1); %after RFT onset
        t2=dsearchn(RFT_time',-0.4); %before cue
        
        for freq=1:2
            scaling_left_low(freq)=mean([cur_RFT.target_left{1,freq}(t1:t2) ; cur_RFT.target_left{3,freq}(t1:t2) ; cur_RFT.dist_left{3,freq}(t1:t2) ; cur_RFT.dist_left{4,freq}(t1:t2)]);
            scaling_right_low(freq)=mean([cur_RFT.target_right{1,freq}(t1:t2) ; cur_RFT.target_right{3,freq}(t1:t2) ;  cur_RFT.dist_right{3,freq}(t1:t2) ;  cur_RFT.dist_right{4,freq}(t1:t2)]);
            scaling_left_high(freq)=mean([cur_RFT.target_left{2,freq}(t1:t2) ; cur_RFT.target_left{4,freq}(t1:t2) ; cur_RFT.dist_left{1,freq}(t1:t2) ; cur_RFT.dist_left{2,freq}(t1:t2)]);
            scaling_right_high(freq)=mean([cur_RFT.target_right{2,freq}(t1:t2) ; cur_RFT.target_right{4,freq}(t1:t2) ;  cur_RFT.dist_right{1,freq}(t1:t2) ;  cur_RFT.dist_right{2,freq}(t1:t2)]);
        end
        
        %scale targets
        for freq=1:2
            for l=1:4
                switch l
                    case {1,3}%low target noise
                        RFT_target_left{l,freq} = (cur_RFT.target_left{l,freq}./scaling_left_low(freq));
                        RFT_target_right{l,freq} = (cur_RFT.target_right{l,freq}./scaling_right_low(freq));
                        RFT_dt_target_left{l,freq} = (cur_RFT_dt.target_left{l,freq}./scaling_left_low(freq));
                        RFT_dt_target_right{l,freq} = (cur_RFT_dt.target_right{l,freq}./scaling_right_low(freq));
                        
                    case {2,4}
                        RFT_target_left{l,freq} = (cur_RFT.target_left{l,freq}./scaling_left_high(freq));
                        RFT_target_right{l,freq} = (cur_RFT.target_right{l,freq}./scaling_right_high(freq));
                        RFT_dt_target_left{l,freq} = (cur_RFT_dt.target_left{l,freq}./scaling_left_high(freq));
                        RFT_dt_target_right{l,freq} = (cur_RFT_dt.target_right{l,freq}./scaling_right_high(freq));
                end
            end
        end
        
        %scale distractors
        for freq=1:2
            for l=1:4
                switch l
                    case {1,2}%high distractor noise
                        RFT_dist_left{l,freq} = (cur_RFT.dist_left{l,freq}./scaling_left_high(freq));
                        RFT_dist_right{l,freq} = (cur_RFT.dist_right{l,freq}./scaling_right_high(freq));
                        RFT_dt_dist_left{l,freq} = (cur_RFT_dt.dist_left{l,freq}./scaling_left_high(freq));
                        RFT_dt_dist_right{l,freq} = (cur_RFT_dt.dist_right{l,freq}./scaling_right_high(freq));
                    case {3,4}%low distractor noise
                        RFT_dist_left{l,freq} = (cur_RFT.dist_left{l,freq}./scaling_left_low(freq));
                        RFT_dist_right{l,freq} = (cur_RFT.dist_right{l,freq}./scaling_right_low(freq));
                        RFT_dt_dist_left{l,freq} = (cur_RFT_dt.dist_left{l,freq}./scaling_left_low(freq));
                        RFT_dt_dist_right{l,freq} = (cur_RFT_dt.dist_right{l,freq}./scaling_right_low(freq));
                end
            end
        end                
        
        %We use a mean relative baseline per noisy level and
        %hemifield, similar to the RFT normalization. This should line up
        %better with the indices as well.
        
        %Get baseline (250ms away from stim onset & Cue onset)
        t1=dsearchn(alpha_time',-1.1);%250 after face onset
        t2=dsearchn(alpha_time',-0.6);%250ms before cue onset
                
        scaling_left_low=mean([cur_alpha.target_left{1}(t1:t2) ; cur_alpha.target_left{3}(t1:t2) ; cur_alpha.dist_left{3}(t1:t2) ; cur_alpha.dist_left{4}(t1:t2)]);
        scaling_right_low=mean([cur_alpha.target_right{1}(t1:t2) ; cur_alpha.target_right{3}(t1:t2) ;  cur_alpha.dist_right{3}(t1:t2) ;  cur_alpha.dist_right{4}(t1:t2)]);
        scaling_left_high=mean([cur_alpha.target_left{2}(t1:t2) ; cur_alpha.target_left{4}(t1:t2) ; cur_alpha.dist_left{1}(t1:t2) ; cur_alpha.dist_left{2}(t1:t2)]);
        scaling_right_high=mean([cur_alpha.target_right{2}(t1:t2) ; cur_alpha.target_right{4}(t1:t2) ;  cur_alpha.dist_right{1}(t1:t2) ;  cur_alpha.dist_right{2}(t1:t2)]);
        
        %scale targets
        for l=1:4
            switch l
                case {1,3}%low target noise
                    Alpha_target_left{l} = (cur_alpha.target_left{l}./scaling_left_low);
                    Alpha_target_right{l} = (cur_alpha.target_left{l}./scaling_right_low);
                    Alpha_dt_target_left{l} = (cur_alpha_dt.target_left{l}./scaling_left_low);
                    Alpha_dt_target_right{l} = (cur_alpha_dt.target_left{l}./scaling_right_low);
                case{2,4} %high target noise
                    Alpha_target_left{l} = (cur_alpha.target_left{l}./scaling_left_high);
                    Alpha_target_right{l} = (cur_alpha.target_left{l}./scaling_right_high);
                    Alpha_dt_target_left{l} = (cur_alpha_dt.target_left{l}./scaling_left_high);
                    Alpha_dt_target_right{l} = (cur_alpha_dt.target_left{l}./scaling_right_high);
            end
        end
        
        %scale distractors
        for l=1:4
            switch l
                case {1,2}%high distractor noise
                    Alpha_dist_left{l} = (cur_alpha.dist_left{l}./scaling_left_high);
                    Alpha_dist_right{l} = (cur_alpha.dist_left{l}./scaling_right_high);
                    Alpha_dt_dist_left{l} = (cur_alpha_dt.dist_left{l}./scaling_left_high);
                    Alpha_dt_dist_right{l} = (cur_alpha_dt.dist_left{l}./scaling_right_high);
                case{3,4} %low distractor noise
                    Alpha_dist_left{l} = (cur_alpha.dist_left{l}./scaling_left_low);
                    Alpha_dist_right{l} = (cur_alpha.dist_left{l}./scaling_right_low);
                    Alpha_dt_dist_left{l} = (cur_alpha_dt.dist_left{l}./scaling_left_low);
                    Alpha_dt_dist_right{l} = (cur_alpha_dt.dist_left{l}./scaling_right_low);
            end
        end
                
        %collapse data
        for l=1:4
            RFT_target{d,l}=mean([RFT_target_left{l,1} RFT_target_left{l,2} RFT_target_right{l,1} RFT_target_right{l,2}],2);
            RFT_dist{d,l}=mean([RFT_dist_left{l,1} RFT_dist_left{l,2} RFT_dist_right{l,1} RFT_dist_right{l,2}],2);
            RFT_dt_target{d,l}=mean([RFT_dt_target_left{l,1} RFT_dt_target_left{l,2} RFT_dt_target_right{l,1} RFT_dt_target_right{l,2}],2);
            RFT_dt_dist{d,l}=mean([RFT_dt_dist_left{l,1} RFT_dt_dist_left{l,2} RFT_dt_dist_right{l,1} RFT_dt_dist_right{l,2}],2);
            
            Alpha_target{d,l}=mean([Alpha_target_left{l} Alpha_target_right{l}],2);
            Alpha_dist{d,l}=mean([Alpha_dist_left{l} Alpha_dist_right{l}],2);
            Alpha_dt_target{d,l}=mean([Alpha_dt_target_left{l} Alpha_dt_target_right{l}],2);
            Alpha_dt_dist{d,l}=mean([Alpha_dt_dist_left{l} Alpha_dt_dist_right{l}],2);
        end
        
        d=d+1;
    else
        disp(['Omitting subject ' int2str(n)])
    end
end

%% Save processed timecourses
timecourses.alpha.ami=ami_alpha;
timecourses.alpha.ami_dt=ami_alpha_dt;
timecourses.alpha.target=Alpha_target;
timecourses.alpha.target_dt=Alpha_dt_target;
timecourses.alpha.dist=Alpha_dist;
timecourses.alpha.dist_dt=Alpha_dt_dist;
timecourses.alpha.LMI_dist_salient=Alpha_dist_LMI_salient_dist;
timecourses.alpha.LMI_dist_noisy=Alpha_dist_LMI_noisy_dist;
timecourses.alpha.LMI_dt_dist_salient=Alpha_dt_dist_LMI_salient_dist;
timecourses.alpha.LMI_dt_dist_noisy=Alpha_dt_dist_LMI_noisy_dist;
timecourses.alpha.DE_low_load=Alpha_target_DE_low_load;
timecourses.alpha.DE_high_load=Alpha_target_DE_high_load;
timecourses.alpha.DE_dt_low_load=Alpha_dt_target_DE_low_load;
timecourses.alpha.DE_dt_high_load=Alpha_dt_target_DE_high_load;
timecourses.alpha.DE_dist_low_load=Alpha_dist_DE_low_load;
timecourses.alpha.DE_dt_dist_low_load=Alpha_dt_dist_DE_low_load;
timecourses.alpha.DE_dist_high_load=Alpha_dist_DE_high_load;
timecourses.alpha.DE_dt_dist_high_load=Alpha_dt_dist_DE_high_load;

timecourses.RFT.target=RFT_target;
timecourses.RFT.target_dt=RFT_dt_target;
timecourses.RFT.dist=RFT_dist;
timecourses.RFT.dist_dt=RFT_dt_dist;
timecourses.RFT.LMI_dist_salient=RFT_dist_LMI_salient_dist;
timecourses.RFT.LMI_dist_noisy=RFT_dist_LMI_noisy_dist;
timecourses.RFT.LMI_dt_dist_salient=RFT_dt_dist_LMI_salient_dist;
timecourses.RFT.LMI_dt_dist_noisy=RFT_dt_dist_LMI_noisy_dist;
timecourses.RFT.DE_low_load=RFT_target_DE_low_load;
timecourses.RFT.DE_high_load=RFT_target_DE_high_load;
timecourses.RFT.DE_dt_low_load=RFT_dt_target_DE_low_load;
timecourses.RFT.DE_dt_high_load=RFT_dt_target_DE_high_load;

timecourses.time=alpha_time;
timecourses.dt_time=alpha_dt_time;

save([group_folder 'Timecourses.mat'],'timecourses');



