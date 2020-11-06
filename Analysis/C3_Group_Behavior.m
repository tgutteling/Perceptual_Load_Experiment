function [] = C3_Group_Behavior()
%
% Agregate behavior at group level
%

%set path
run('/rds/projects/2017/jenseno-01/Tjerk/set_path');

%folders
proc_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/';

%Check for High frequency TFRs in subject folders (these also contain the
%behavioral data)
sub_folders=dir([proc_folder filesep 'S*']);
cnt=1;
for s=1:size(sub_folders,1)
    if exist([proc_folder sub_folders(s).name filesep sub_folders(s).name '_all_clean.mat'])>0
        datasets{cnt}=[proc_folder sub_folders(s).name filesep sub_folders(s).name '_all_clean.mat'];
        cnt=cnt+1;
    end
end

%load data
disp(['Found ' int2str(length(datasets)) ' processed datasets'])
for d=1:length(datasets)
    disp(['Loading ' datasets{d}])
    cur=load(datasets{d});
    cur=cur.data;
    % behavior structure:
    % cfg.behavior(i,1) = i;           %trailNr
    % cfg.behavior(i,2) = ConMat(i,1); %left or right attention
    % cfg.behavior(i,3) = ConMat(i,2); %condition
    % cfg.behavior(i,4) = ConMat(i,4); %left or right eye movement
    % cfg.behavior(i,5) = RespKey;     %left or right key response
    % cfg.behavior(i,6) = (figureframes-cueframes)*ifi;     %delay time in s
    % cfg.behavior(i,7) = rt;          %in ms
    % cfg.behavior(i,8) = rt2;          %in ms
    % cfg.behavior(i,9) = rt3;          %in ms
    % cfg.behavior(i,10) = correct;     %answer
    
    % conditions:
    %1 - low target load, low distractor load
    %2 - high target load, low distractor load
    %3 - low target load, high distractor load
    %4 - high target load, high distractor load
    
    % RT per condition:
    for c=1:2
        for l=1:4                  
            beh{d}.RT.left(c,l)=nanmean(cur{c}.left{l}.behavior(cur{c}.left{l}.behavior(:,10)==1,8))*1000;
            beh{d}.RT.right(c,l)=nanmean(cur{c}.right{l}.behavior(cur{c}.right{l}.behavior(:,10)==1,8))*1000;
            beh{d}.RT.all(c,l)=nanmean([cur{c}.left{l}.behavior(cur{c}.left{l}.behavior(:,10)==1,8) ; cur{c}.right{l}.behavior(cur{c}.right{l}.behavior(:,10)==1,8)])*1000;
            beh{d}.RT.raw{c,l}=[cur{c}.left{l}.behavior ; cur{c}.right{l}.behavior];
            
            beh{d}.hits.left(c,l)=nansum(cur{c}.left{l}.behavior(:,10))/length(cur{c}.left{l}.behavior(:,10))*100;
            beh{d}.hits.right(c,l)=nansum(cur{c}.right{l}.behavior(:,10))/length(cur{c}.right{l}.behavior(:,10))*100;
            beh{d}.hits.all(c,l)=nansum([cur{c}.left{l}.behavior(:,10) ; cur{c}.right{l}.behavior(:,10)])/length([cur{c}.left{l}.behavior(:,10) ; cur{c}.right{l}.behavior(:,10)])*100;
            
            beh{d}.miss.left(c,l)=sum(isnan(cur{c}.left{l}.behavior(:,10)))/length(cur{c}.left{l}.behavior(:,10))*100;
            beh{d}.miss.right(c,l)=sum(isnan(cur{c}.right{l}.behavior(:,10)))/length(cur{c}.right{l}.behavior(:,10))*100;
            beh{d}.miss.all(c,l)=sum(isnan([cur{c}.left{l}.behavior(:,10) ; cur{c}.right{l}.behavior(:,10)]))/length([cur{c}.left{l}.behavior(:,10) ; cur{c}.right{l}.behavior(:,10)])*100;                        
        end
    end
    
end

save([proc_folder 'group' filesep 'group_Behavior.mat'],'beh');

end