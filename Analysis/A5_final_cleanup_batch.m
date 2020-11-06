function [] = A5_final_cleanup_batch(subject,av_type)
%
% Post ICA visual inspection of remaining artifacts
% After cleanup, data are sorted in conditions and saved


%% Set up folders

%set path
run('/rds/projects/2017/jenseno-01/Tjerk/set_path');

%processed data location (for saving)
data_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/';
layout_folder='/rds/projects/2017/jenseno-01/Tjerk/Load2/Layouts/';

switch av_type
    case 1 %Cue-locked
        if exist([data_folder filesep subject filesep subject '_icclean.mat'])==0
            error('No ICA cleaned data available for this subject!')
        end
        
        fprintf('Loading data..')
        load([data_folder filesep subject filesep subject '_icclean.mat']);
        disp('Done.')
        
    case 2 %Discrimination target-locked
        if exist([data_folder filesep subject filesep subject '_icclean_dt.mat'])==0
            error('No ICA cleaned data available for this subject!')
        end
        
        fprintf('Loading data..')
        load([data_folder filesep subject filesep subject '_icclean_dt.mat']);
        disp('Done.')
end

sub_nr=str2num(subject(2:end));

%% First step: Trial rejection
%Results of visual inspection is noted down in the subjects_rejection(/_dt)
%mat files, specifying trials to be rejected. 
switch av_type
    case 1
        load([data_folder 'subject_rejection.mat']);
    case 2
        load([data_folder 'subject_rejection_dt.mat']);
end

%% Second step: use eyelink to identify eye movement trials

switch av_type
    case 1
        %Subject should fixate centrally at least between t=[-1 -0.3]
        t1=dsearchn(data.el.time',-1);
        t2=dsearchn(data.el.time',-0.3);
        el=data.el;
    case 2
        %There is no fixation period, but subject should fixate before DT onset
        t1=dsearchn(data.el_dt.time',-2);
        t2=dsearchn(data.el_dt.time',-1.1);
        el=data.el_dt;
end
x_center=median(el.x(:,t1:t2)');
y_center=median(el.y(:,t1:t2)');

%Make [0 0] the center
el_x=(el.x-repmat(x_center',1,size(el.x,2)));
el_y=(el.y-repmat(y_center',1,size(el.y,2)));

%We have a screen distance and size, we can calculate degrees/pixel
width=data.exp_cfg{1}.width;                 %projection screen width in cm
height=data.exp_cfg{1}.height;                %projection screen height in cm
dist=data.exp_cfg{1}.dist;

deg_p=atand((width/1920)/dist);

%set threshold in pixels
thres_deg=3; %it's about 3degrees from fixation to the edge of the stimulus
time_thres=25; %minimum time duration of eye movement to be rejected (avoids glitches)
pix_thres=thres_deg/deg_p;

%now let's check trials with above threshold eye movements
switch av_type
    case 1
        t1=dsearchn(el.time',-1);
        t2=dsearchn(el.time',1);
    case 2
        t1=dsearchn(el.time',-1.5);
        t2=dsearchn(el.time',0);
end

el_bad=unique(sort([find((sum(abs(el_x(:,t1:t2))>pix_thres,2))>time_thres) ; find((sum(abs(el_y(:,t1:t2))>pix_thres,2))>time_thres)]))';

disp([int2str(length(el_bad)) ' trials exceed eye movement threshold'])
disp(['Eye movement trials: ' int2str(el_bad)])

%% Third step: use EOG to identify unfortunate blinks
if av_type==1
    EOG_chan=find(strcmp('vEOG',data.meg.label));
    time=data.meg.time{1};
else
    EOG_chan=find(strcmp('vEOG',data.meg_dt.label));
    time=data.meg_dt.time{1};
end
cfg=[];
cfg.channel=EOG_chan;
if av_type==1
    EOG=ft_selectdata(cfg,data.meg);
else
    EOG=ft_selectdata(cfg,data.meg_dt);
end
cfg=[];
cfg.keeptrials='yes';
EOG=ft_timelockanalysis(cfg,EOG);
EOG=squeeze(EOG.trial);

%so, S01 has impedance check on trial 1:57
if strcmp(subject,'S01')
    EOG(1:57,:)=NaN;
end

%Correct scaling for session 1 and 2
%get the cutoff point
end_trl=find(data.config==data.config(1),1,'last'); %last trial of first block 

%get resting EOG
med1=nanmean(median(EOG(1:end_trl,:),2));
med2=nanmean(median(EOG(end_trl+1:end,:),2));

%find peak/blink values
peaks1=max(abs(EOG(1:end_trl,:)),[],2);peaks1(isnan(peaks1))=[];
sorted_peaks1=sort(peaks1);
midway1=((mean(sorted_peaks1(end-5:end))-mean(sorted_peaks1(1:25)))/2);
max1=median(peaks1(peaks1>midway1));
peaks2=max(abs(EOG(end_trl+1:end,:)),[],2);peaks2(isnan(peaks2))=[];
sorted_peaks2=sort(peaks2);
midway2=((mean(sorted_peaks2(end-5:end))-mean(sorted_peaks2(1:25)))/2);
max2=median(peaks2(peaks2>midway2));

%remove mean and apply peak scaling
EOG_z=NaN(size(EOG));
EOG_z(1:end_trl,:)=(EOG(1:end_trl,:)-med1)./max1;
EOG_z(end_trl+1:end,:)=(EOG(end_trl+1:end,:)-med2)./max2;

%EOG_thres_z=3;
EOG_thres_z=0.75;

%now let's check trials with above threshold blinks
if av_type==1
    t1=dsearchn(time',-1);
    t2=dsearchn(time',1);
else
    t1=dsearchn(time',-1.5);
    t2=dsearchn(time',0.25);
end

%z-transform
eog_bad=find(sum(abs(EOG_z(:,t1:t2))>EOG_thres_z,2))';

disp([int2str(length(eog_bad)) ' trials exceed blink threshold'])
disp(['Blink trials: ' int2str(eog_bad)])


%% Remove bad trials from data
cfg=[];
bad=unique(sort([reject_major_tech_art{sub_nr} reject_major_sub_art{sub_nr} el_bad eog_bad]));
trials=1:length(data.config);
trials(bad)=[];
cfg.trials=trials;
if av_type==1
    data.meg=ft_selectdata(cfg,data.meg);
else
    data.meg_dt=ft_selectdata(cfg,data.meg_dt);
end

disp([int2str(length(bad)) ' trials rejected, ' int2str(length(trials)) ' trials left. ' num2str((length(trials)/512)*100,3) '% of the data preserved.'])

%and additional data
data.config(bad)=[];
data.behavior(bad,:)=[];
data.ConMat(bad,:)=[];

%eyelink
switch av_type
    case 1
        data.el.config(bad)=[];
        data.el.att_side(bad)=[];
        data.el.BlockNr(bad)=[];
        data.el.TrialNr(bad)=[];
        data.el.x(bad,:)=[];
        data.el.y(bad,:)=[];
        data.el.pupil(bad,:)=[];
    case 2
        data.el_dt.config(bad)=[];
        data.el_dt.att_side(bad)=[];
        data.el_dt.BlockNr(bad)=[];
        data.el_dt.TrialNr(bad)=[];
        data.el_dt.x(bad,:)=[];
        data.el_dt.y(bad,:)=[];
        data.el_dt.pupil(bad,:)=[];
end

%% Channel interpolation

%Get neighbours
load([layout_folder 'Neighbours_grad2.mat']);
load([layout_folder 'Neighbours_grad3.mat']);
load([layout_folder 'Neighbours_mag.mat']);

%channel selection
switch av_type
    case 1
        MEG_sens=strmatch('MEG',data.meg.label);
        sens_type=str2num(cellfun(@(x) x(end),data.meg.label(MEG_sens),'UniformOutput',1));
        nsens=length(data.meg.label);
    case 2
        MEG_sens=strmatch('MEG',data.meg_dt.label);
        sens_type=str2num(cellfun(@(x) x(end),data.meg_dt.label(MEG_sens),'UniformOutput',1));
        nsens=length(data.meg_dt.label);
end

mags=sort([MEG_sens(sens_type==1)]);
grad2_sens=sort([MEG_sens(sens_type==2)]);
grad3_sens=sort([MEG_sens(sens_type==3)]);
aux=setdiff(1:nsens,MEG_sens);

%split up datasets
switch av_type
    case 1
        cfg=[];
        cfg.channel=mags;
        tmp_mags=ft_selectdata(cfg,data.meg);tmp_mags.grad=data.grad_common;%data.grad{1};%
        cfg.channel=grad2_sens;
        tmp_grad2=ft_selectdata(cfg,data.meg);tmp_grad2.grad=data.grad_common;%data.grad{1};%
        cfg.channel=grad3_sens;
        tmp_grad3=ft_selectdata(cfg,data.meg);tmp_grad3.grad=data.grad_common;%data.grad{1};%
        cfg.channel=aux;
        tmp_aux=ft_selectdata(cfg,data.meg);
    case 2
        cfg=[];
        cfg.channel=mags;
        tmp_mags=ft_selectdata(cfg,data.meg_dt);tmp_mags.grad=data.grad_common;%data.grad{1};
        cfg.channel=grad2_sens;
        tmp_grad2=ft_selectdata(cfg,data.meg_dt);tmp_grad2.grad=data.grad_common;%data.grad{1};
        cfg.channel=grad3_sens;
        tmp_grad3=ft_selectdata(cfg,data.meg_dt);tmp_grad3.grad=data.grad_common;%data.grad{1};
        cfg.channel=aux;
        tmp_aux=ft_selectdata(cfg,data.meg_dt);
end

%channel interpolation per sensor type and direction
cfg=[];
cfg.method='spline';

%mag interpolation
mag_missing=setdiff({neighbours_mag.label}',tmp_mags.label);
if ~isempty(mag_missing)
    cfg.neighbours=neighbours_mag;
    cfg.missingchannel=mag_missing;
    tmp_mags=ft_channelrepair(cfg,tmp_mags);
end

%grad2 interpolation
grad2_missing=setdiff({neighbours_grad2.label}',tmp_grad2.label);
if ~isempty(grad2_missing)
    cfg.neighbours=neighbours_grad2;
    cfg.missingchannel=grad2_missing;
    tmp_grad2=ft_channelrepair(cfg,tmp_grad2);
end

%grad3 interpolation
grad3_missing=setdiff({neighbours_grad3.label}',tmp_grad3.label);
if ~isempty(grad3_missing)
    cfg.neighbours=neighbours_grad3;
    cfg.missingchannel=grad3_missing;
    tmp_grad3=ft_channelrepair(cfg,tmp_grad3);    
end

if av_type==1
    data.meg=ft_appenddata([],tmp_mags,tmp_grad2,tmp_grad3,tmp_aux);
    data.meg.grad=data.grad_common;
else
    data.meg_dt=ft_appenddata([],tmp_mags,tmp_grad2,tmp_grad3,tmp_aux);
    data.meg_dt.grad=data.grad_common;
end


%% Sorting
%Now that we have all data, cleaned up nicely, lets sort conditions
if av_type==1
    if sum(data.ConMat(:,2)~=data.el.config')>0
        error('Mismatch between MEG data and eyelink file!')
    end
else
    if sum(data.ConMat(:,2)~=data.el_dt.config')>0
        error('Mismatch between MEG data and eyelink file!')
    end
end

%first sort by configuration
order=unique(data.config,'stable');
cfg=[];
switch  av_type
    case 1
        for i=1:2
            cfg.trials=find(data.config==i);
            c{i}.meg=ft_selectdata(cfg,data.meg);
            c{i}.behavior=data.behavior(cfg.trials,:);
            c{i}.ConMat=data.ConMat(cfg.trials,:);
            
            %eyelink
            c{i}.el.x=data.el.x(cfg.trials,:);
            c{i}.el.y=data.el.y(cfg.trials,:);
            c{i}.el.pupil=data.el.pupil(cfg.trials,:);
            c{i}.el.time=data.el.time;
            
            %additional info
            c{i}.grad=data.grad{order(i)};
            c{i}.exp_cfg=data.exp_cfg{order(i)};
        end
        
        %now separate left and right attention trials and sort per LOAD condition
        tmp=[];
        for config=1:2
            for i=1:4
                cfg=[];

                %Attention left
                cfg.trials=find(c{config}.ConMat(:,1)==1 & c{config}.ConMat(:,2)==i);
                tmp{config}.left{i}.meg=ft_selectdata(cfg,c{config}.meg);
                tmp{config}.left{i}.behavior=c{config}.behavior(cfg.trials,:);
                tmp{config}.left{i}.ConMat=c{config}.ConMat(cfg.trials,:);
                tmp{config}.left{i}.el.x=c{config}.el.x(cfg.trials,:);
                tmp{config}.left{i}.el.y=c{config}.el.y(cfg.trials,:);
                tmp{config}.left{i}.el.pupil=c{config}.el.pupil(cfg.trials,:);

                %Attention right
                cfg.trials=find(c{config}.ConMat(:,1)==2 & c{config}.ConMat(:,2)==i);
                tmp{config}.right{i}.meg=ft_selectdata(cfg,c{config}.meg);
                tmp{config}.right{i}.behavior=c{config}.behavior(cfg.trials,:);
                tmp{config}.right{i}.ConMat=c{config}.ConMat(cfg.trials,:);
                tmp{config}.right{i}.el.x=c{config}.el.x(cfg.trials,:);
                tmp{config}.right{i}.el.y=c{config}.el.y(cfg.trials,:);
                tmp{config}.right{i}.el.pupil=c{config}.el.pupil(cfg.trials,:);
            end
            tmp{config}.grad=c{config}.grad;
            tmp{config}.exp_cfg=c{config}.exp_cfg;
        end
        
    case 2
        for i=1:2
            cfg.trials=find(data.config==i);
            c{i}.meg_dt=ft_selectdata(cfg,data.meg_dt);
            c{i}.behavior=data.behavior(cfg.trials,:);
            c{i}.ConMat=data.ConMat(cfg.trials,:);
            
            %eyelink
            c{i}.el_dt.x=data.el_dt.x(cfg.trials,:);
            c{i}.el_dt.y=data.el_dt.y(cfg.trials,:);
            c{i}.el_dt.pupil=data.el_dt.pupil(cfg.trials,:);
            c{i}.el_dt.time=data.el_dt.time;
            
            %additional info
            c{i}.grad=data.grad{order(i)};
            c{i}.exp_cfg=data.exp_cfg{order(i)};
        end
        
        %now separate left and right attention trials and sort per LOAD condition
        tmp=[];
        for config=1:2
            for i=1:4
                cfg=[];

                %Attention left
                cfg.trials=find(c{config}.ConMat(:,1)==1 & c{config}.ConMat(:,2)==i);
                tmp{config}.left{i}.meg_dt=ft_selectdata(cfg,c{config}.meg_dt);
                tmp{config}.left{i}.behavior=c{config}.behavior(cfg.trials,:);
                tmp{config}.left{i}.ConMat=c{config}.ConMat(cfg.trials,:);
                tmp{config}.left{i}.el_dt.x=c{config}.el_dt.x(cfg.trials,:);
                tmp{config}.left{i}.el_dt.y=c{config}.el_dt.y(cfg.trials,:);
                tmp{config}.left{i}.el_dt.pupil=c{config}.el_dt.pupil(cfg.trials,:);

                %Attention right
                cfg.trials=find(c{config}.ConMat(:,1)==2 & c{config}.ConMat(:,2)==i);
                tmp{config}.right{i}.meg_dt=ft_selectdata(cfg,c{config}.meg_dt);
                tmp{config}.right{i}.behavior=c{config}.behavior(cfg.trials,:);
                tmp{config}.right{i}.ConMat=c{config}.ConMat(cfg.trials,:);
                tmp{config}.right{i}.el_dt.x=c{config}.el_dt.x(cfg.trials,:);
                tmp{config}.right{i}.el_dt.y=c{config}.el_dt.y(cfg.trials,:);
                tmp{config}.right{i}.el_dt.pupil=c{config}.el_dt.pupil(cfg.trials,:);
            end
            tmp{config}.grad=c{config}.grad;
            tmp{config}.exp_cfg=c{config}.exp_cfg;
        end        
end

%replace old data structure with new
data = [];
data=tmp;

%% Saving
fprintf(['Saving data for ' subject '...'])
if ~exist([data_folder subject],'dir')>0    
    disp('Creating subject directory..')
    mkdir([data_folder subject]);    
end
if av_type==1
    save([data_folder subject filesep subject '_all_clean.mat'],'data','-v7.3');
else
    save([data_folder subject filesep subject '_all_clean_dt.mat'],'data','-v7.3');
end
    
fprintf('..Done\n')

end
