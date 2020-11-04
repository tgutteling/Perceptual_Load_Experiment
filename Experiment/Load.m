function [] = Load()

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Experiment with Rapid Tagging
% 18/10/2018
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load Experiment

% A pop up screen appears when running the script,
% Initials: first block normal, second block add 2
% if block = 2, the freq tagging configuration is swapped

% Press e during the experiment to stop eyetracker feedback
% Press q during the experiment to quit

% 22/08/2018 - added practice mode
% 23/08/2018 - added intro
% 17/09/2018 - corrected luminance for all stimuli to 50%
% 24/09/2018 - added nata button box (untested)
% 26/09/2018 - added 'shuffle' noise mode to maintain luminance levels across loads
% 28/09/2018 - second implementation of NAtA button boxes using KbQueue
% 16/10/2018 - replaced KbqueueWait with Check because it did not seem to work as well
% 18/10/2018 - added option to use button boxes in practice mode
% 18/10/2018 - Modded eyelink check to make it easier
% 04/11/2020 - Cleaned up script and added comments

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization

% Clear the workspace and the screen
close all;
clear all;
sca;
fclose('all');
Screen('Preference', 'SkipSyncTests', 1); %must be 0 during experiment
AssertOpenGL;
KbName('UnifyKeyNames'); % for easy use of Keyboard keys
PsychDefaultSetup(2);    % call some default settings for setting up Psychtoolbox

cfg = [];

%Physical screen parameters (cm) (defaults)
cfg.width=70.6;                 %projection screen width in cm
cfg.height=49.8;                %projection screen height in cm
cfg.dist=147.5;                    %distance from subject eye to screen in cm

%add images folder to path (assuming script is in the same folder as images folder)
p = mfilename('fullpath');
exp_dir=p(1:find(p==filesep,1,'last'));
addpath([p(1:find(p==filesep,1,'last')) 'Images'])
cfg.exp_dir=exp_dir;

prompt = {'Subject code:', ...
    'Block Number:', ...
    'Practice Mode:', ...
    'Frequency Counterbalance: ', ...
    'Screen width (cm): ', ...
    'Screen height (cm): ', ...
    'Screen distance (cm): '};

dlg_title = 'Input';
num_lines = 1;
defaultans = {'Test','1','0','1',num2str(cfg.width),num2str(cfg.height),num2str(cfg.dist)};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
cfg.sub=answer{1};
cfg.block=str2num(answer{2});

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN settings

%Here is the specification of all parameters used in the study

%screen setup
cfg.PracticeMode = str2num(answer{3});%Enabling practice mode means single screen mode (non-rapid), no flicker, half contrast, no eyelink
cfg.RunIntro=0;                   %run intro by default?
cfg.rapidMode = 1;                %turn on quadRGB mode for 12 simultaneous frames on the ProPixx
cfg.use_screenres = 0;            %use the detected screen resolution
cfg.fullscreen = 0;               %open a fullscreen window, if not, manual resolution will be used
cfg.debugmode = 1;                %for rapidmode testing without triggers
cfg.DataPixxOnly = 0;             %for rapidmode testing without triggers with propixx
cfg.manual_resx = 1920/2;         %manual x resolution
cfg.manual_resy =1080/2;          %manual y resolution
cfg.TextSize=12;                  %Default text size, will scale to other screen resolutions
cfg.photoDiode=1;                 %Add corner stimulus for photodiode?
cfg.diodeSize=2;                  %Size of the photodiode in cm

%stimuli
cfg.StimSize=8;                  %Stimulus size in vis ang
cfg.FixationSize = 10;           %Size of the central fixation point, in Pixels
cfg.FixationsShift=20;           %Upward shift of screen center, to allow more space below fixation, in percent (0-100)
cfg.eccentricity=7;              %Eccentricity of stimuli, in deg vis ang

%Physical screen parameters (cm) (input)
cfg.width=str2num(answer{5});   %projection screen width in cm
cfg.height=str2num(answer{6});  %projection screen height in cm
cfg.dist=str2num(answer{7});    %distance from subject eye to screen in cm

%Frequencies
f1=63;                          %Frequencies to be used in frequency tagging, will be counterbalanced
f2=70;
cfg.WaveShape='sin';            %Shape of the waveform 'sin' for sinusoidal, 'square' for square wave
cfg.Phaselock=1;                %1=Phase locked RFT, same phase RFT stimulation for every trial; 0=random phase
cfg.FreqBins=8;                 %Nr of frequency bins with non-phase locked RFT
cfg.diodeFreq=3;                %Frequency to present to photodiode? 1=f1, 2=f2, 3=left stimulus, depending on config

%Timing
cfg.base_t = 1;                 %Duration (s) of the baseline period grey background
cfg.cue_t = 0.35;               %Duration (s) of the cueing period
cfg.stim_t = 1;                 %Duration (s) of the Figure stimulation (oriented backgroud+figure)
cfg.stim_jitter = 1;            %Baseline jitter. Stimulus duration is stim_t+jitter*[0..1]
cfg.move_t = 0.15;              %Found in yeshurun paper 2013
cfg.maxresp = 1;                %Max duration (s) of the response period
cfg.feedback1 = 1;              %Duration (s) of the feedback
cfg.iti = 1;                    %intertrial interval
cfg.BreakTrials=26;             %How many trials between breaks?

%stimuli & block
cfg.noiseType='shuffle';         %Type of noise. 'shuffle'=randomize pixel positions, keeps luminance equal;'salt': Salt & pepper noise, adds black and white pixels
cfg.HighNoiseDensity = 0.5;     %Target noise density, this affects approximately d*numel(I) pixels. Behavior was tested using 0.35 Salt&Pepper noise
cfg.LowNoiseDensity = 0;         %Distr noise density, zero density to create non noised images
cfg.cond_reps = 1;               %Nr of repetitions of unique conditions per block

%Eyelink
cfg.el.eyelink=0;               %eyelink on/off
cfg.el.Eyeused='LEFT_EYE';      %eye used
cfg.el.edffile = [cfg.sub(1:3) '_' int2str(cfg.block)];  %EDF filename, taken from subject data and session nr
cfg.el.fixation_window=4;       %Size of the fixation window for online fixation control (deg vis ang) - Self--> 2 deg
cfg.el.feedback=1;              %Feedback (arrows) on (1) or off (0);
cfg.el.eyelinkKey=KbName('e');  %Key used to toggle eyelink feedback on or off during experiment.

%initialize keys
cfg.KeyLeft=KbName('4$'); %should key NAtA button box left index (5 button version)
cfg.KeyRight=KbName('7&'); %should key NAtA button box right index (5 button version)
cfg.quitKey=KbName('q');
cfg.continueKey=KbName('c');
cfg.escapeKey = KbName('ESCAPE');

%block preparation
%set initial frequency
cfg.CounterBalance=str2num(answer{4});
switch str2num(answer{4}) %counterbalance over subjects
    case 1
        cfg.FreqMat = [f1 f2];        
    case 2
        cfg.FreqMat = [f2 f1];
end

%set block (session) frequency counterbalance over blocks
if ~mod(str2num(answer{2}),2) %even block, i.e. block 2, 4 etc
    cfg.FreqMat=circshift(cfg.FreqMat,[0 1]); %swap left/right frequencies
end

if str2num(answer{2})>1
    cfg.RunIntro=0;
end

%Set parameters for practice mode
if cfg.PracticeMode
    cfg.rapidMode = 0;
    cfg.use_screenres = 1;
    cfg.fullscreen = 1;
    cfg.debugmode = 1;
    cfg.photoDiode=0;
    cfg.FreqMat= [0 0];
    cfg.el.eyelink=0;
    cfg.Phaselock=1;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LogFile
%Initialize a text logfile to be written to during experiment. 

% use log fid for logging experiment procedures and keeping track were errors appear
timestamp=clock;
stampstr=[int2str(timestamp(1)) '-' int2str(timestamp(2)) '-' int2str(timestamp(3)) '-' int2str(timestamp(4)) '-' int2str(timestamp(5)) '-' int2str(timestamp(6))];

%Create logging directory
if ~exist([exp_dir 'Logs' filesep cfg.sub],'dir')>0
    mkdir([exp_dir 'Logs' filesep cfg.sub]);
end

%Create logfile
fid=fopen([exp_dir 'Logs' filesep cfg.sub filesep 'Log_' cfg.sub '_ses_' int2str(cfg.block) '_' stampstr '.txt'],'w');
cfg.fid=fid;
WaitSecs(1);
fprintf(fid,['This logfile was created on ' datestr(clock) '\n']);

%Log the configuration options
fprintf(fid,'Config options: \n');
names=fieldnames(cfg);
vals=struct2cell(cfg);
for i=1:length(names)
    if ~isstruct(vals{i})
        if ~isstr(vals{i})
            fprintf(fid,[names{i} ': ' num2str(vals{i}) '\n']);
        else
            fprintf(fid,'%s',[names{i} ': ' vals{i}]);
            fprintf(fid,'\n');
        end
    else
        names2=fieldnames(getfield(cfg,names{i}));
        vals2=struct2cell(getfield(cfg,names{i}));
        for i2=1:length(names2)
            fprintf(fid,[names2{i2} ': ' num2str(vals2{i2}) '\n']);
        end
    end
end

cfg.LogFile=[cfg.sub filesep 'Log_' cfg.sub '_ses_' int2str(cfg.block) '_' stampstr '.txt'];
cfg.datestamp=stampstr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Screen setup

% Get the screen numbers
screens =Screen('Screens');

%select screen
screenNumber = 1;%max(screens);
cfg.screenNumber = screenNumber;
log(fid,['Using screen ' int2str(screenNumber)]);

% Get the size of the on screen window and set resolution
sc_info=Screen('Resolution', screenNumber);
screenXpixels=sc_info.width;
screenYpixels=sc_info.height;

log(fid,['Measured screen size is ' int2str(screenXpixels) 'x' int2str(screenYpixels)]);
if cfg.use_screenres
    log(fid,'Using measured screen resolution');
    resx=screenXpixels;
    resy=screenYpixels;
else
    log(fid,['Using manual screen resolution of ' int2str(cfg.manual_resx) 'x' int2str(cfg.manual_resy)]);
    resx=cfg.manual_resx;
    resy=cfg.manual_resy;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Eyelink
%Set up the eyelink for the experiment
el_override=0; %No eyelink in actual experiment, use only in case of fault

if cfg.el.eyelink
    %add eyelink script folder (should be in main experiment folder)
    addpath([exp_dir filesep 'Eyelink']);
    
    %make directory if it doesn't already exist (local computer)
    cfg.el.eyedir = [exp_dir filesep 'Eyelink' filesep ];
    if ~exist(cfg.el.eyedir, 'dir'); mkdir(cfg.el.eyedir);end
    
    %check whether files already exist for this subject/session
    if exist([exp_dir 'Eyelink' filesep cfg.el.edffile '.edf'])>0
        cont=input('Warning! Eyelink file will be overwritten, do you want to continue? (y/n) ','s');
        if cont=='n'
            error('Session aborted')
        end
    end
    
    cfg.el_rect = [0 0 resx resy];
    
    % Set parameters, start and calibrate eyelink
    cfg = el_Start(cfg);
    
    %Parameters for fixation control, note conversion to pixels does not
    %take into account rapidmode, because eyelink knows only projected screen
    fixWinSize=round(((tand(cfg.el.fixation_window)*cfg.dist)*(resy/cfg.height)));
    cfg.el.fixationWindow = [-fixWinSize -fixWinSize fixWinSize fixWinSize];
    
else
    if ~cfg.debugmode %is this is real experiment time, eyelink should be on
        if el_override
            warning('Eyelink not in use, continuing anyway...')
            log(fid,'Eyelink NOT initialized, default overriden');
        else
            cont=input('Warning, Eyelink is not selected! do you want to continue without eyelink? (y/n) ','s');
            if cont=='n'
                error('Experiment aborted')
            elseif cont=='y'
                warning('Eyelink not in use, continuing anyway...')
                log(fid,'Eyelink NOT selected, default overriden');
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parallel Port IO & triggers %

% Send triggers for:
% - trial start w type
% - fig onset w location
% - discrimination target onset
% - manual response

% set up triggers
TriggerStart = 1;   % Trial start
TriggerFig1 = 2;    % config 1  / Low target load, low distractor salience
TriggerFig2 = 4;    % config 2  / High target load, low distractor salience
TriggerFig3 = 8;    % config 3  / Low target load, high distractor salience
TriggerFig4 = 16;   % config 4  / High target load, high distractor salience
TriggerMov = 32;    % Onset of stimulus eye movement
TriggerResp = 64;   % subject response
triggers=[TriggerFig1 TriggerFig2 TriggerFig3 TriggerFig4];

log(fid,['Trigger ' int2str(TriggerStart) ' is Trial start']);
log(fid,['Trigger ' int2str(TriggerFig1) ' is Config Low Low']);
log(fid,['Trigger ' int2str(TriggerFig2) ' is Config High Low']);
log(fid,['Trigger ' int2str(TriggerFig3) ' is Config Low High']);
log(fid,['Trigger ' int2str(TriggerFig4) ' is Config High High']);
log(fid,['Trigger ' int2str(TriggerResp) ' is subject Response']);

if ~cfg.debugmode    
    PortAddress = hex2dec('BFF8');
    ioObjTrig = io64;
    status = io64(ioObjTrig);
    io64(ioObjTrig,PortAddress,0); %trigger 0 (reset)
    
    cfg.PortAddress=PortAddress;
    cfg.ioObjTrig=ioObjTrig;
else
    log(fid,'DebugMode, not sending triggers!')
end

cfg.Triggers=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Open screen

% Open an on-screen window if not fullscreen
if  cfg.fullscreen
    log(fid,'Opening full screen window');
    [window, ~] = PsychImaging('OpenWindow', screenNumber, 0.5); % 0.5 makes gray
else
    log(fid,['Opening a ' int2str(cfg.manual_resx) 'x' int2str(cfg.manual_resy) ' window']);
    offset=50;
    [window, ~] = PsychImaging('OpenWindow', screenNumber, 0.5, [offset offset cfg.manual_resx+offset cfg.manual_resy+offset]);
end
cfg.window=window;
%enable alpha blending
log(fid,'Enabling alpha blending');
Screen(window,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%Hide the cursus when fullscreen
if  cfg.fullscreen
    HideCursor;
end

%Flip to clear
Screen('Flip', window);

%Query the frame duration
ifi = Screen('GetFlipInterval', window);
log(fid,['Measured interflip interval: ' num2str(ifi,5) '. Measured refresh rate: ' num2str(1/ifi,5)]);

%Set maximum priority level
topPriorityLevel = MaxPriority(window);
if ~cfg.debugmode
    Priority(topPriorityLevel);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stimulus matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ConMat stimulus matrix structure
%column 1: direction, attend left or right (cue direction)
% 1 = left
% 2 = right
%column 2: configuration, (low and high load target and distractor)
% 1 = low load target, low salience distractor
% 2 = high load target, low salience distractor
% 3 = low load target, high salience distractor
% 4 = high load target, high salience distractor
%column 3: faces (eight different faces)
%column 4: eye movements, left or right switch (1 = left, 2 = right)
%column 5: distractor movements: 1 = congruent 2 = incongruent.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%total number of unique conditions
direc = 2;      % Left and Right Attention
config = 4;     % high and low target and distractor loads
face = 8;       % number of different faces
conds = direc*config*face; % total amount of conditions
ConMat = NaN(conds,3);     % create NaN matrix for all conditions to be filled

log(fid,['Generating condition matrix with ' int2str(conds*4) ' unique conditions']);

%Create matrix with all combinations
direc_cnt = 1;
config_cnt = 1;
face_cnt = 2;
dir_chk = conds/direc;
conf_chk = dir_chk/config;
ConMat(1,:) = [1 1 1];
for i = 2:conds
    if mod(i-1,conf_chk) == 0
        config_cnt = config_cnt+1;
        face_cnt = 1;
    end
    if mod(i-1,dir_chk) == 0
        direc_cnt = direc_cnt+1;
        config_cnt = 1;
        face_cnt = 1;
    end
    ConMat(i,:) = [direc_cnt,config_cnt,face_cnt];
    face_cnt = face_cnt+1;
end

%Replicate matrix, to add extra left and right eye movements condition
ConMat = repmat(ConMat,2,1);
ConMat(:, 4) = [repmat(1, conds, 1); repmat(2,conds, 1)];

%Replicate matrix, to add extra congruent vs incongruent conditions
%In which 1 is congruent and 2 is incongruent
ConMat=repmat(ConMat,2,1);
ConMat(:, 5) = [repmat(1, 2*conds, 1); repmat(2,2*conds, 1)];

%For multiple repetitions of unique conditions, replicate matrix (now 1)
ConMat=repmat(ConMat,cfg.cond_reps,1);

%Randomly permute
ConMat = ConMat(randperm(size(ConMat,1)),:);

%save in cfg
cfg.ConMat = ConMat;

%total nr of trials per block
nTrials = size(ConMat,1);
log(fid,['There will be ' int2str(nTrials) ' trials per block']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Stimulus Placement
% Due to the propixx rapidmode, stimulus placement is more complicated, as
% it has to be replicated over 4 quadrants per screen
log(fid,'Calculating stimulus positions');

%Define the size of the stimuli, in pixels, derived from visual angle
if cfg.rapidMode
    stim_size=round((tand(cfg.StimSize)*cfg.dist)*(resy/(cfg.height*2)));
else
    stim_size=round(((tand(cfg.StimSize)*cfg.dist)*(resy/cfg.height)));
end
log(fid,['Using a ' int2str(stim_size) ' pixel width for the stimuli']);

%Calculate the stimulus centers
xPos=resx/4;
yPos=resy/4;
pos_1 = round([xPos yPos]);
pos_2 = round([3*(xPos) yPos]);
pos_3 = round([xPos 3*(yPos)]);
pos_4 = round([3*(xPos) 3*(yPos)]);
qcenters =[pos_1 ; pos_2 ; pos_3 ; pos_4]; % centers for each quandrant
q_upshift_pix=round(yPos*(cfg.FixationsShift*1e-2));
nm_upshift_pix=round((resy/2)*(cfg.FixationsShift*1e-2)); %shift of quadrant center in pixels

if cfg.rapidMode
    log(fid,['Shifting the central fixation upward by ' int2str(cfg.FixationsShift) 'percent, or ' int2str(q_upshift_pix) ' pixels'])
else
    log(fid,['Shifting the central fixation upward by ' int2str(cfg.FixationsShift) 'percent, or ' int2str(nm_upshift_pix) ' pixels'])
end

%apply screen center shift
qcenters(:,2)=qcenters(:,2)-q_upshift_pix;
center_nm=[resx/2 (resy/2)-nm_upshift_pix];
cfg.center_nm=center_nm;
cfg.qcenters=qcenters;

%for every quandrant, calculate rects (useful for text placement)
q_rects(1,:)=[0 0 qcenters(1,:)*2];
q_rects(2,:)=[resx/2 0 resx qcenters(2,2)*2];
q_rects(3,:)=[0 resy/2 resx/2 resy-q_upshift_pix*2];
q_rects(4,:)=[resx/2 resy/2 resx resy-q_upshift_pix*2];
cfg.q_rects=q_rects;

%Get the placement of the stimuli
ecc_xy=sqrt(cfg.eccentricity^2/2); %calculate x and y shift, given eccentricity, euclidean (eccentricity is measured as distance from fixation point to stimulus center)
if cfg.rapidMode
    ecc = round(((tand(ecc_xy)*cfg.dist)*((resy/2)/cfg.height)));
    for q=1:4
        positions_stim{q}=[qcenters(q,:)+[-ecc ecc] ; qcenters(q,:)+[ecc ecc]];
    end
    if cfg.photoDiode
        
        %calculate size of photodiode in pixels
        diode_size_pix=round(0.5*(resy/cfg.height)*cfg.diodeSize);
        
        %calculate diode positions (bottom right)
        diode_pos{1}=[resx/2-diode_size_pix resy/2-diode_size_pix resx/2 resy/2];
        diode_pos{2}=[resx-diode_size_pix resy/2-diode_size_pix resx resy/2];
        diode_pos{3}=[resx/2-diode_size_pix resy-diode_size_pix resx/2 resy];
        diode_pos{4}=[resx-diode_size_pix resy-diode_size_pix resx resy];
    end
else
    %NORMAL MODE
    ecc = round(((tand(ecc_xy)*cfg.dist)*(resy/cfg.height))); %euclidean
    positions_stim=[(resx/2)-ecc (resy/2)+ecc-nm_upshift_pix ; (resx/2)+ecc (resy/2)+ecc-nm_upshift_pix];
    if cfg.photoDiode
        %calculate size of photodiode in pixels
        diode_size_pix=round((resy/cfg.height)*cfg.diodeSize);
        
        %calculate diode positions (bottom right)
        diode_pos=[resx-diode_size_pix resy-diode_size_pix resx resy];
    end
end

% load feedback figure
load([exp_dir filesep 'Images' filesep 'arrow4_sm.mat']); %fixation reminder arrow

% Placement of ARROW for Fixation
s=size(arrow_image);
rect=[0 0 s(1)/(2*1920/resx) s(2)/(2*1080/resy)]; %scales relative to 'standard' 1920x1080
FixArrPos=cell(1,5);
FixArrPos_nm = floor(CenterRectOnPointd(rect,resx/2,resy/2-nm_upshift_pix)); %'normal mode' position for feedback arrow

arr_rect=[0 0 s(1)/(4*1920/resx) s(2)/(4*1080/resy)];
for q=1:4
    FixArrPos{q} = floor(CenterRectOnPointd(arr_rect,qcenters(q,1),qcenters(q,2)));
end
arrow = Screen('MakeTexture', window, arrow_image);

%Placement of eyelink fixation window, note that it's in the actual center
%of (projected) screen, even in rapidmode
if cfg.el.eyelink
    cfg.el.fixationWindow = floor(CenterRectOnPointd(cfg.el.fixationWindow,resx/2,resy/2-nm_upshift_pix));
end

%Set the text size, default is good for 1920x1080 in rapidmode (i.e. 540 pixels in y direction)
if cfg.rapidMode
    cfg.TextSize=cfg.TextSize*resy/1080;
else
    cfg.TextSize=cfg.TextSize*resy/540;
end
Screen('TextSize', window, cfg.TextSize);

%FXIATION CROSS and Cueing lines
%Set size of the arms and linewidth
if cfg.rapidMode
    if cfg.fullscreen
        fixCrossDimPix = 10; lineWidthPix = 2;
    else
        fixCrossDimPix = 5; lineWidthPix = 1;
    end
else %larger cues in normal mode
    if cfg.fullscreen
        fixCrossDimPix = 20; lineWidthPix = 4;
    else
        fixCrossDimPix = 10; lineWidthPix = 2;
    end
end

cfg.lineWidthPix=lineWidthPix;

%Get coordinates
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
CrossCoords = [xCoords; yCoords];
cfg.CrossCoords=CrossCoords;

xCoords = [-fixCrossDimPix fixCrossDimPix -fixCrossDimPix fixCrossDimPix];
yCoords = [fixCrossDimPix 0 -fixCrossDimPix 0];
RightCueCoords = [xCoords; yCoords];

xCoords = [-fixCrossDimPix fixCrossDimPix -fixCrossDimPix fixCrossDimPix];
yCoords = [0 -fixCrossDimPix 0 fixCrossDimPix];
LeftCueCoords = [xCoords; yCoords];

cfg.RightCueCoords=RightCueCoords;
cfg.LeftCueCoords=LeftCueCoords;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure Stimuli

log(fid,'Creating loaded and non-loaded face stimuli');

% Textures:
% faceTex(facenr,variation: LL(center, left, right), HL(center, left, right)

% load mask figure
mask = rgb2gray(imread('Mask.jpg')); % same size as stimuli, 128 for surround and 0 for center

% create alpha layer with the size of the stimulus, to make the
% placeholder (rest of square around the face) transparant
% Putting the alpha layer as the second layer of the image will do so
imageSize = size(mask);
ci = [199,199,197];
[xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
mask_alpha = uint8((xx.^2 + yy.^2)<ci(3)^2);
mask_alpha(mask_alpha==1)=255;  % Alpha layer: center is 255, surrounding is 0

if cfg.photoDiode
    diode=uint8(ones(size(mask_alpha))*255);
    diode(:,:,2)=mask_alpha;
    diode_tex = Screen('MakeTexture', window, diode);
end

faceNr = [2 3 4 5 6 7 8 9]; %for loading the right faces from the Images folder

% faceTex matrix:
% rows: face variations
% columns:
%     1 - neutral face        high noise
%     2 - left gazing face    high noise
%     3 - right gazing face   high noise
%     4 - neutral face        low noise
%     5 - left gazing face    low noise
%     6 - right gazing face   low noise
log(fid,['Using noise level of ' num2str(cfg.LowNoiseDensity) ' for Low-noise faces'])
log(fid,['Using noise level of ' num2str(cfg.HighNoiseDensity) ' for High-noise faces'])

% loop over faces to create faceTex matrix
for i = 1:length(faceNr)
    
    f = faceNr(i);  %get the right faceNr.
    
    %For eacht stimulus loop over high and low load to create textures for both conditions
    for loads = 1 : 2 %1=high, 2=low noise
        if loads == 1
            dens = cfg.HighNoiseDensity; %use prespecified noise density
            jit = 0;                 %save loaded images on first three columns
        else
            dens = cfg.LowNoiseDensity;
            jit = 3;                 %save non-loaded images on 3rd 4th and 5th matrix column
        end
        
        % load images with CENTER configuration
        centerconfig = sprintf('f%i_center2.png', f);       %create name variable
        img = imread(centerconfig);
        img = rgb2gray(img);
        
        if loads==1 && strcmp(cfg.noiseType,'shuffle')
            %add noise by shuffling pixels
            %img = rgb2gray(img);
            FaceOnly=img(mask_alpha==255);
            randPixels=round(dens*size(FaceOnly,1));
            PixelSelect = (round(rand(randPixels,1)*(length(FaceOnly)-1)))+1;
            FaceOnly(PixelSelect)=Shuffle(FaceOnly(PixelSelect));
            img(mask_alpha==255)=FaceOnly;
        else if strcmp(cfg.noiseType,'salt') && loads==1
                %add salt and pepper noise and create grey scaled image
                %Mask the image to calculate the luminance of the face only
                img = imnoise(img,'salt & pepper',dens);
                masked_img=(mask_alpha/255).*img;
                %Equalize luminance
                if mean(masked_img(:))/255>0.5 %too much white, reduce lights
                    img=((0.5/(mean(masked_img(:))/255)))*img;
                else %too dark, increase darks
                    img=255-((mean(masked_img(:))/255)/0.5)*(255-img);
                end
            end
        end
                
        % add alpha mask on the second layer
        img(:,:,2) = mask_alpha;        
        faceTex(i,1+jit) = Screen('MakeTexture', window, img); %save loaded texture in matrix
        img = flip(img,2);
        FlipfaceTex(i,1+jit)  = Screen('MakeTexture', window, img); %save loaded texture in matrix
        
        % load images with LEFT configuration
        leftconfig = sprintf('f%i_left2.png', f);
        img = imread(leftconfig);
        img = rgb2gray(img);
        
        if loads==1 && strcmp(cfg.noiseType,'shuffle')
            %add noise by shuffling pixels
            FaceOnly=img(mask_alpha==255);
            randPixels=round(dens*size(FaceOnly,1));
            PixelSelect = (round(rand(randPixels,1)*(length(FaceOnly)-1)))+1;
            FaceOnly(PixelSelect)=Shuffle(FaceOnly(PixelSelect));
            img(mask_alpha==255)=FaceOnly;
        else if strcmp(cfg.noiseType,'salt') && loads==1
                %add salt and pepper noise and create grey scaled image
                %Mask the image to calculate the luminance of the face only
                img = imnoise(img,'salt & pepper',dens);
                masked_img=(mask_alpha/255).*img;
                %Equalize luminance
                if mean(masked_img(:))/255>0.5 %too much white, reduce lights
                    img=((0.5/(mean(masked_img(:))/255)))*img;
                else %too dark, increase darks
                    img=255-((mean(masked_img(:))/255)/0.5)*(255-img);
                end
            end
        end
        
        img(:,:,2) = mask_alpha;
        faceTex(i, 2+jit) = Screen('MakeTexture', window, img); %left
        img = flip(img,2);
        FlipfaceTex(i,2+jit)  = Screen('MakeTexture', window, img); %save loaded texture in matrix
        
        % load images with RIGHT configuration
        rightconfig = sprintf('f%i_right2.png',f);
        img = imread(rightconfig);
        img = rgb2gray(img);
        if loads==1 && strcmp(cfg.noiseType,'shuffle')
            %add noise by shuffling pixels
            FaceOnly=img(mask_alpha==255);
            randPixels=round(dens*size(FaceOnly,1));
            PixelSelect = (round(rand(randPixels,1)*(length(FaceOnly)-1)))+1;
            FaceOnly(PixelSelect)=Shuffle(FaceOnly(PixelSelect));
            img(mask_alpha==255)=FaceOnly;
        else if strcmp(cfg.noiseType,'salt') && loads==1
                %add salt and pepper noise and create grey scaled image
                %Mask the image to calculate the luminance of the face only                
                img = imnoise(img,'salt & pepper',dens);
                masked_img=(mask_alpha/255).*img;
                %Equalize luminance
                if mean(masked_img(:))/255>0.5 %too much white, reduce lights
                    img=((0.5/(mean(masked_img(:))/255)))*img;
                else %too dark, increase darks
                    img=255-((mean(masked_img(:))/255)/0.5)*(255-img);
                end
            end
        end
        
        img(:,:,2) = mask_alpha;
        faceTex(i, 3+jit) = Screen('MakeTexture', window, img); %right
        img = flip(img,2);
        FlipfaceTex(i,3+jit)  = Screen('MakeTexture', window, img); %save loaded texture in matrix
    end
end

%store in cfg (for intro)
cfg.faceTex=faceTex;

%now that we have the stimuli, we can calculate the placement in PTB terms (rects)
patchrect=[0 0 stim_size stim_size];%
if cfg.rapidMode
    positions=cell(1,4);
    for q=1:4
        for p=1:2
            positions{q}(p,:)=floor(CenterRectOnPointd(patchrect,positions_stim{q}(p,1),positions_stim{q}(p,2)));
        end
    end
else
    positions=zeros(2,4);
    for p=1:2
        positions(p,:)=floor(CenterRectOnPointd(patchrect,positions_stim(p,1),positions_stim(p,2)));
    end
end
cfg.positions=positions;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial KeyBoard settings
KbName('UnifyKeyNames');

if ~cfg.debugmode
    log(fid,['Using button box ' KbName(cfg.KeyLeft) ' (left) and ' KbName(cfg.KeyRight) ' (right) as response buttons.'])
else
    cfg.BBLeft=cfg.KeyLeft;
    cfg.BBRight=cfg.KeyRight;
    cfg.KeyLeft = KbName('LeftArrow');
    cfg.KeyRight = KbName('RightArrow');
end

%set active keys (others will be ignored)
if cfg.PracticeMode
    cfg.ActiveKeys=[cfg.BBLeft cfg.BBRight KbName('LeftArrow') KbName('RightArrow') cfg.quitKey cfg.el.eyelinkKey cfg.escapeKey cfg.continueKey];
else
    cfg.ActiveKeys=[cfg.KeyLeft cfg.KeyRight cfg.quitKey cfg.el.eyelinkKey cfg.escapeKey cfg.continueKey];
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Frequency Tagging Timecourse
%The Rapid frequency tagging timecourse is precomputed according to the
%maximal amount of frames in a trial.
% calculate the amount of frames needed for each part of the experiment.

% small recall of durations from basic settings
% Timing settings:
% cfg.base_t          %Duration (s) of the baseline period
% cfg.cue_t           %Duration (s) of the cueing period
% cfg.stim_t          %Duration (s) of the Figure stimulation
% cfg.stim_jitter     %Baseline jitter. Stimulus duration is stim_t+jitter*[0..1]
% cfg.move_t          %Found in yeshurun paper 2013
% cfg.maxresp         %Max duration (s) of the response period
% cfg.feedback1       %Duration (s) of the feedback

%Get the maximal amount of frames to calculate the timecourse for
base_f = round(cfg.base_t/ifi);
cue_f = round(cfg.cue_t/ifi);
stim_f = round(cfg.stim_t/ifi);
jitter_frames = round(cfg.stim_jitter/ifi);
move_f = round(cfg.move_t/ifi);
maxresp_f = round(cfg.maxresp/ifi);

feedback1_f = round(cfg.feedback1/ifi);
max_trialframes = round((cfg.base_t+cfg.cue_t+cfg.stim_t+cfg.stim_jitter+cfg.maxresp)/ifi);%max frames per trial

%Log amounts of frames
log(fid,'Calculating presentation frames');
log(fid,['Using ' int2str(base_f) ' baseline frames']);
log(fid,['Using ' int2str(cue_f) ' cueing frames']);
log(fid,['Using ' int2str(stim_f) ' to ' int2str(stim_f+jitter_frames) ' stimulus frames']);
log(fid,['Using maximally ' int2str(maxresp_f) ' frames for participants response']);
log(fid,['Using ' int2str(feedback1_f) ' feedback frames if trial is missed']);
log(fid,['Maximal amount of trial frames: ' int2str(max_trialframes)]);

%adjust times to an integer number of frames, to avoid rounding differences
%difference will be marginal
cfg.base_t = base_f*ifi;
cfg.cue_t = cue_f*ifi;
cfg.stim_t = stim_f*ifi;
cfg.stim_jitter = jitter_frames*ifi;
cfg.maxresp = maxresp_f*ifi;

%Log times per state
log(fid,'Based on the measured presentation rate, actual times will be:');
log(fid,['Baseline: ' num2str(cfg.base_t,7) 's']);
log(fid,['Cue: ' num2str(cfg.cue_t,7) 's']);
log(fid,['Stimulation: ' num2str(cfg.stim_t,7) 's']);
log(fid,['Jitter: ' num2str(cfg.stim_jitter,7) 's']);
log(fid,['Max response time:' num2str(cfg.maxresp,7) 's']);
log(fid,['Max trial length:' num2str(max_trialframes*ifi,7) 's']);

%Here we define the phase timecourse of the frequencies. To enable trial
%averaging in the time domain, phase will be 0 at t=0 (zeroframe), stimulus (figure)
%onset
if cfg.rapidMode
    frame_mult=12; %every refresh is 12 frames
else
    frame_mult=1; %or just one
end

%Effective presentation frequency. Should be 1440 for propixx
Fs=(1/ifi)*frame_mult;
cfg.Fs=Fs;
log(fid,['Effective refresh rate: ' num2str(Fs,6) 'Hz']);

if sum(cfg.FreqMat>Fs/2)>0
    warning('Presentation rate too low for the chosen flicker frequencies!')
end

%Frequency timecourse parameters
cfg.patch_amplitude = 0.5;
cfg.patch_startPhase = 0;
cfg.f_offset = 0;
log(fid,'Frequency parameters: ');
log(fid,['RFT amplitude: ' num2str(cfg.patch_amplitude,2) 's']);
log(fid,['RFT start phase (t=0): ' num2str(cfg.patch_startPhase,2)]);
log(fid,['RFT offset: ' num2str(cfg.f_offset,2)]);

%initialize the table
if cfg.Phaselock
    freqTable=NaN(length(cfg.FreqMat),(max_trialframes*frame_mult));
else
    freqTable=NaN(length(cfg.FreqMat),cfg.FreqBins,(max_trialframes*frame_mult));
end
frametime=NaN(1,frame_mult*max_trialframes);

for f = 1:length(cfg.FreqMat)
    patch_frequency = cfg.FreqMat(f);
    patch_angFreq = 2 * pi * patch_frequency;
    start_time=(cfg.base_t+cfg.cue_t)*-1;
    frametime=start_time:ifi/frame_mult:(max_trialframes*frame_mult)*(ifi/frame_mult)+start_time;
    frametime=frametime(1:max_trialframes*frame_mult);
    if cfg.Phaselock
        if strcmpi(cfg.WaveShape,'square') %square wave
            freqTable(f,:)= cfg.patch_amplitude * square(patch_angFreq * frametime + cfg.patch_startPhase) + cfg.patch_amplitude + cfg.f_offset;
        else %sinusoidal
            freqTable(f,:)= cfg.patch_amplitude * sin(patch_angFreq * frametime + cfg.patch_startPhase) + cfg.patch_amplitude + cfg.f_offset;
        end
    else
        for b=1:cfg.FreqBins
            cfg.patch_startPhase =(2*pi/cfg.FreqBins)*(b-1);
            if strcmpi(cfg.WaveShape,'square') %square wave
                freqTable(f,b,:)= cfg.patch_amplitude * square(patch_angFreq * frametime + cfg.patch_startPhase) + cfg.patch_amplitude + cfg.f_offset;
            else %sinusoidal
                freqTable(f,b,:)= cfg.patch_amplitude * sin(patch_angFreq * frametime + cfg.patch_startPhase) + cfg.patch_amplitude + cfg.f_offset;
            end
        end
    end
end

cfg.freqTable=freqTable;

if cfg.photoDiode
    if cfg.Phaselock
        switch cfg.diodeFreq
            case 1
                cfg.diodeTable=freqTable(dsearchn(cfg.FreqMat',f1),:);
            case 2
                cfg.diodeTable=freqTable(dsearchn(cfg.FreqMat',f2),:);
            case 3
                cfg.diodeTable=freqTable(1,:);
        end
    end
end

%calculate all permutation of phase differences. It is important to balance
%the phase differences well to properly cancel out phase interference
tmp=[];
combs=[];
intervals=[1:cfg.FreqBins]-1;
for i=1:cfg.FreqBins
    tmp(:,1)=[1:cfg.FreqBins]';
    tmp(:,2)=mod(intervals+(i-1),cfg.FreqBins)+1;
    combs=[combs ; tmp];
end
cfg.combs=combs;

if ~cfg.Phaselock
    %randomize phase table for trials
    reps=nTrials/length(combs);
    if mod(nTrials/(length(combs)),1)>0
        warning(['Nr of frequency bins (' int2str(cfg.FreqBins) ') does not fit into an integer nr of trials (' int2str(nTrials) ')'])
    end
    combsMat=repmat(combs,reps,1);
    cfg.combsMat=combsMat(randperm(size(combsMat,1)),:);
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define different states of the EXPERIMENT
%create a specification of states of a trial
%use the exact timing within the syncframes loop to switch between
%states.

log(fid,'Defining states');
trial_states = struct;

%define baseline
trial_states.baseline.name='baseline';
for s=1:2 %two stimuli
    trial_states.baseline.stimuli(s).type = 'Texture';
    trial_states.baseline.stimuli(s).flicker =s ; %0 for no flicker, >0 for position in freqTable
end
% add fixation cross
trial_states.baseline.stimuli(3).type = 'DrawLines';
trial_states.baseline.stimuli(3).coords = {CrossCoords, CrossCoords};
trial_states.baseline.stimuli(3).params = {lineWidthPix, center_nm};
trial_states.baseline.stimuli(3).flicker = 0;

%define 'cue' display
trial_states.cue.name = 'cue';
for s=1:2
    trial_states.cue.stimuli(s).type = 'Texture';
    trial_states.cue.stimuli(s).flicker =s;
end
% add cue left or right
trial_states.cue.stimuli(3).type = 'DrawLines';
trial_states.cue.stimuli(3).coords = {LeftCueCoords, RightCueCoords};
trial_states.cue.stimuli(3).params = {lineWidthPix, center_nm};
trial_states.cue.stimuli(3).flicker = 0;


%define 'delay' display
trial_states.stim_onset.name = 'stim_onset';
for s=1:2
    trial_states.stim_onset.stimuli(s).type = 'Texture';
    trial_states.stim_onset.stimuli(s).flicker = s;
end
% add fixation cross
trial_states.stim_onset.stimuli(3).type = 'DrawLines';
trial_states.stim_onset.stimuli(3).coords = {CrossCoords, CrossCoords};
trial_states.stim_onset.stimuli(3).params = {lineWidthPix, center_nm};
trial_states.stim_onset.stimuli(3).flicker = 0;

%define 'eye movement' display
trial_states.movement.name = 'movement';
for s=1:2
    trial_states.movement.stimuli(s).type = 'Texture';
    trial_states.movement.stimuli(s).flicker = s;
end
% add fixation cross
trial_states.movement.stimuli(3).type = 'DrawLines';
trial_states.movement.stimuli(3).coords = {CrossCoords, CrossCoords};
trial_states.movement.stimuli(3).params = {lineWidthPix, center_nm};
trial_states.movement.stimuli(3).flicker = 0;

%define 'response' display
trial_states.response.name = 'response';
%only a fixation cross
trial_states.response.stimuli(1).type = 'DrawLines';
trial_states.response.stimuli(1).coords = {CrossCoords, CrossCoords};
trial_states.response.stimuli(1).params = {lineWidthPix, center_nm};
trial_states.response.stimuli(1).flicker = 0;

%define 'feedback1' display
trial_states.feedback1.name = 'feedback1';
trial_states.feedback1.stimuli(1).type = 'Text';
trial_states.feedback1.stimuli(1).text = 'Too late, please try to respond a bit faster!';


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXPERIMENTAL LOOP

% loop over trials
% loop over syncframes. (the amount of frames within a trial )

% create states to put information from current trial part in, used to put
% into the screen function later
% switch between the states available

log(fid,'EXPERIMENT START');

% Setup Propixx 1440 Hz
if cfg.rapidMode && ~cfg.debugmode || cfg.DataPixxOnly
    Datapixx('Open');
    Datapixx('SetPropixxDlpSequenceProgram', 5); % 2 for 480, 5 for 1440 Hz, 0 for normal
    Datapixx('RegWrRd');
end

%Run intro?
if cfg.RunIntro
    run_intro(cfg);
end

states={'baseline','cue', 'stim_onset', 'movement', 'response', 'feedback1'};

if cfg.el.eyelink
    %Experiment start message to eyelink
    Eyelink('Message', 'Experiment start');
end

log(fid,['BLOCK: ' int2str(cfg.block)]);

%set frequencies per quadrant - since we're not using multiple blocks now,
BlockFreq = cfg.FreqMat;
BlockFreqTable = freqTable;
log(fid,'Frequencies: ');
log(fid, ['Left: ' int2str(BlockFreq(1))]);
log(fid, ['Right: ' int2str(BlockFreq(2))]);

if cfg.el.eyelink
    %send block start message to eyelink
    Eyelink('Message', ['Block ' int2str(cfg.block) ' start']);
end

%no trigger sent yet
cfg.triggerSent=0;

%Initizalize the KbQueue for listening to very short responses (NAtA button boxes)
keylist=zeros(1,256); %Set all keys to zero (ignore)
keylist(cfg.ActiveKeys)=1; %set active keys to 1 (listen)
KbQueueCreate(0,keylist);%%Create queue


for i = 1 : nTrials
    
    %Store current trial
    cfg.cur_trial=i;
    
    %Start listening to input
    KbQueueStart;
    KbQueueFlush;% clear all keyboard presses so far
       
    %welcome screen at first trial
    if i == 1
        
        if cfg.rapidMode % in rapid mode at different location
            for q = 1:4                
                [~,~,~]=DrawFormattedText(window, 'The experiment is about to start.\nPress any button to START!', 'center', 'center',1,[],[],[],[],[],q_rects(q,:));
            end
        else % normal mode
            DrawFormattedText(window, 'The experiment is about to start.\nPress any key to START!', 'center', center_nm(2),1)
        end
        vbl = Screen('Flip', window);
        GetResp=1;
        while GetResp %turns at eye movement
            [keyIsDown, ~]=KbQueueCheck();
            if keyIsDown
                GetResp=0;
            end
        end
        
        %add empty screen
        Screen('Flip', window);
        WaitSecs(cfg.iti);
        
    else
        %add empty screen
        Screen('Flip', window);
        WaitSecs(cfg.iti);
    end
    
    %Set flags
    GetResp = 0;
    breakFlag = 0;
    rt = NaN;
    rt2 = NaN;
    rt3 = NaN;
    correct = NaN;
    RespKey = NaN;
    
    %Set trial dependent face stimulus
    faceStim = ConMat(i,3); % get the right faceNr from the condition matrix
    
    % get load configurations sorted out for this trial
    
    % faceTex matrix:
    % rows: face variations
    % columns:
    %     1 - neutral face        high noise
    %     2 - left gazing face    high noise
    %     3 - right gazing face   high noise
    %     4 - neutral face        low noise
    %     5 - left gazing face    low noise
    %     6 - right gazing face   low noise
    
    % attend left
    if ConMat(i,1) == 1
        if ConMat(i,2) == 1         % low load target & low distractor salience
            faceTr = [faceTex(faceStim,4) faceTex(faceStim,1)];
            FlipfaceTr = [FlipfaceTex(faceStim,4) FlipfaceTex(faceStim,1)];
        
        elseif ConMat(i,2) == 2     % high load target & low distractor salience
            faceTr = [faceTex(faceStim,1) faceTex(faceStim,1)];
            FlipfaceTr = [FlipfaceTex(faceStim,1) FlipfaceTex(faceStim,1)];
            
        elseif ConMat(i,2) == 3     % low load target & high distractor salience
            faceTr = [faceTex(faceStim,4) faceTex(faceStim,4)];
            FlipfaceTr = [FlipfaceTex(faceStim,4) FlipfaceTex(faceStim,4)];
            
        elseif ConMat(i,2) == 4     % high load target & high distractor salience
            faceTr = [faceTex(faceStim,1) faceTex(faceStim,4)];
            FlipfaceTr = [FlipfaceTex(faceStim,1) FlipfaceTex(faceStim,4)];            
        end
    else % attend right
        if ConMat(i,2) == 1         % low target & low distractor
            faceTr = [faceTex(faceStim,1) faceTex(faceStim,4)];
            FlipfaceTr = [FlipfaceTex(faceStim,1) FlipfaceTex(faceStim,4)];
        
        elseif ConMat(i,2) == 2     % high low
            faceTr = [faceTex(faceStim,1) faceTex(faceStim,1)];
            FlipfaceTr = [FlipfaceTex(faceStim,1) FlipfaceTex(faceStim,1)];
            
        elseif ConMat(i,2) == 3     % low high
            faceTr = [faceTex(faceStim,4) faceTex(faceStim,4)];
            FlipfaceTr = [FlipfaceTex(faceStim,4) FlipfaceTex(faceStim,4)];
            
        elseif ConMat(i,2) == 4     % high high
            faceTr = [faceTex(faceStim,4) faceTex(faceStim,1)];
            FlipfaceTr = [FlipfaceTex(faceStim,4) FlipfaceTex(faceStim,1)];
            
        end
    end
    
    NewfaceTR = [FlipfaceTr(1) faceTr(2)];
    % for left and right faces use the same variable but add 1 or 2 to
    % get the same loaded face with eyes left or right
    
    %define frames for each state and transition
    baseframes = base_f;
    cueframes = baseframes + cue_f;
    figureframes = cueframes + stim_f + round((cfg.stim_jitter*rand)/ifi); %determine amount of jitter.
    movement = figureframes + move_f;
    response = movement + 1; %maxresp_f; %no frames are drawn during response
    feedback1 = response + feedback1_f;
    syncframes = feedback1;
    
    %save exact delay time for each trial.
    figureTime = figureframes*ifi; % in s
    
    % create list containing all states, we will loop over this using state_count
    state_transitions=[baseframes cueframes figureframes movement response feedback1];
    
    %set frequency table
    if cfg.Phaselock
        cur_freqTable=BlockFreqTable;
    else %randomize phase
        comb_sel=cfg.combsMat(i,:);
        
        bin1=comb_sel(1);
        bin2=comb_sel(2);
        cur_freqTable(1,:)=squeeze(BlockFreqTable(1,bin1,:))';
        cur_freqTable(2,:)=squeeze(BlockFreqTable(2,bin2,:))';
        
        switch cfg.diodeFreq
            case 1
                cfg.diodeTable=cur_freqTable(dsearchn(cfg.FreqMat',f1),:);
            case 2
                cfg.diodeTable=cur_freqTable(dsearchn(cfg.FreqMat',f2),:);
            case 3
                cfg.diodeTable=cur_freqTable(1,:);
        end
    end
    
    %lets start with baseline
    cur_state = trial_states.baseline;
    state_trig = TriggerStart; %trial start (baseline)
    state_change=1; %we just changed the state, we want a flip timestamp
    
    %Trigger to be used in this trial, indicating condition
    StimTrig = triggers(ConMat(i,2));
    
    %Define eye movement face
    %check for incongruent or congruent trial
    if ConMat(i,5) == 2 % incongruent
        if ConMat(i,1) == 1 %attend left
            if ConMat(i,4) == 1 % eyemove left
                % links links rechts rechts
                MovefaceTR = [FlipfaceTr(1)+4 faceTr(2)+4];
            else % == 2  %eyemove right
                %rechts links
                MovefaceTR = [FlipfaceTr(1)+2 faceTr(2)+2];
            end
        elseif ConMat(i,1) == 2 %attend right
            if ConMat(i,4) == 1 %eyemove left
                % rechts links
                MovefaceTR = [FlipfaceTr(1)+2 faceTr(2)+2];
            else % == 2  %rechts
                % links rechts
                MovefaceTR = [FlipfaceTr(1)+4 faceTr(2)+4];
            end
        end
        
    else        % CONGRUENT!: eye movement towards same side
        if ConMat(i,1) == 1 || ConMat(i,1) ==2 %attend left en right
            if ConMat(i,4) == 1 %links
                % left left
                MovefaceTR = [FlipfaceTr(1)+4 faceTr(2)+2];
            else % == 2  %right
                % right right
                MovefaceTR = [FlipfaceTr(1)+2 faceTr(2)+4];
            end
        end
    end
    
    %With EYELINK enabled, before trial start, subject should be fixating
    fix = 0;
    blink = 0;
    fix_loss = 0;
    arr_on = 0;
    el_toggle=0;
    KbQueueFlush;
    if cfg.el.eyelink && cfg.el.feedback
        %Draw standard background and fixation while waiting
        if cfg.rapidMode
            for q=1:4 %four quadrants draw fixation cross
                Screen('DrawLines', window, CrossCoords, lineWidthPix, 1, qcenters(q,:), 2);
            end
        else %normal mode
            Screen('DrawLines', window, CrossCoords, lineWidthPix, 1, center_nm, 2);
        end
        % Flip to the screen
        vbl = Screen('Flip', window);
        
        while fix==0
            %check whether fixation is ok
            sample = Eyelink('NewestFloatSample');
            % Get current gaze position from sample
            x = sample.gx(1); %first sample should be left eye
            y = sample.gy(1);
            if x~=cfg.el.defaults.MISSING_DATA && y~=cfg.el.defaults.MISSING_DATA && sample.pa(1)>0 %check if its an actual measurement
                %compare x and y to fixation window
                fix = x > cfg.el.fixationWindow(1) &&  x <  cfg.el.fixationWindow(3) && y > cfg.el.fixationWindow(2) && y < cfg.el.fixationWindow(4);
            end
            
            %Here we check whether a key is pressed that will abort the
            %experiment (to interrupt without crashing matlab)            
            [ keyIsDown, ~, keyCode ] = KbCheck;
            
            % Check pressed key, if its the escape key (e.g. 'q'), abort the
            % experiment.
            if keyIsDown && keyCode(cfg.escapeKey)
                cleanup(cfg);
                break;
            end
            
            if keyIsDown && keyCode(cfg.el.eyelinkKey)>0
                if cfg.el.feedback
                    cfg.el.feedback=0;
                    log(fid,'Disabling eyelink feedback');
                    break;
                else
                    cfg.el.feedback=1;
                    log(fid,'Enabling eyelink feedback');
                end
            end
        end
    end
    
    %With that out of the way, let's start
    
    %trial start timestamp
    state_count=1;          %counter used to loop over stimulus states
    t_trial = GetSecs;      %trial start time
    cfg.t_trial=t_trial;
    
    %get rid of lingering keypresses
    KbQueueFlush;
    
    %trial start
    if cfg.el.eyelink
        %send trial start trigger to eyelink
        Eyelink('Message', [int2str(cfg.block) ' ' int2str(cfg.cur_trial) ' Config: ' int2str(ConMat(i,2)) ' Attention side: ' int2str(ConMat(i,1)) ' FaceNr: ' int2str(ConMat(i,3))  ' Movement direction: ' int2str(ConMat(i,4)) ' Congruency: ' int2str(ConMat(i,5))]);
    end
    
    %log info
    log(fid, [int2str(cfg.block) ' ' int2str(cfg.cur_trial) ' Config: ' int2str(ConMat(i,2)) ' Attention side: ' int2str(ConMat(i,1)) ' FaceNr: ' int2str(ConMat(i,3))  ' Movement direction: ' int2str(ConMat(i,4)) ' Congruency: ' int2str(ConMat(i,5))]);
    
    for j = 1:syncframes
        
        %CHANGE STATE
        %transition if amount of frames is achieved & Set appropriate parameters
        if j > state_transitions(state_count) && strcmp(cur_state.name, states{state_count})
            state_count = state_count+1;
            cur_state = getfield(trial_states,states{state_count});
            state_change = 1; % for flip timestamp
            state_trig = StimTrig;
            
            %moving into response state
            if strcmp(cur_state.name, 'response')
                GetResp = 1;
            end
            
            %moving into eye movement state
            if strcmp(cur_state.name, 'movement')
                GetResp = 1;
                NewfaceTR = MovefaceTR;
            end
            
            %Move to feedback state
            if strcmp(cur_state.name, 'feedback1')
                GetResp = 0;
            end
        end
        
        %CHECK FIXATION
        %We'll check eye position during the experiment, once per screen refresh (e.g. 120Hz)
        if cfg.el.eyelink %after initial period
            
            %Get eyelink sample
            sample = Eyelink('NewestFloatSample');
            
            % Get current gaze position from sample
            x = sample.gx(1); %first sample should be left eye
            y = sample.gy(1);
            fix_prev=fix;
            blink_prev=blink;
            if x~=cfg.el.defaults.MISSING_DATA && y~=cfg.el.defaults.MISSING_DATA && sample.pa(1)>0 %check if its an actual measurement
                %compare x and y to fixation window
                fix = x > cfg.el.fixationWindow(1) &&  x <  cfg.el.fixationWindow(3) && y > cfg.el.fixationWindow(2) && y < cfg.el.fixationWindow(4);
                blink=0;
            else
                fix = 1; %eye closed, but not necessarily violating fixation
                blink=1;
            end
            if fix_prev==0 && ~fix && ~fix_loss
                log(fid,[int2str(i) ' FIXATION LOST'])
                start_frame=j;
                fix_loss=1;
            end
            if blink && ~blink_prev
                log(fid,[int2str(i) ' FIXATION LOST'])
            end
            if fix
                fix_loss=0;
                arr_on=0;
            end
        end
        
        %DRAW STIMULI
        if cfg.rapidMode
            for q = 1:4 %for all quadrants
                for p = 1:length(cur_state.stimuli) %for all stimuli
                    % Draw fix cross or cue
                    if strcmpi(cur_state.stimuli(p).type,'DrawLines')
                        Screen('DrawLines', window, cur_state.stimuli(p).coords{ConMat(i,1)}, lineWidthPix, 1, qcenters(q,:), 2);
                    end
                    % Draw feedback text
                    if strcmpi(cur_state.stimuli(p).type,'Text')
                        [~,~,~]=DrawFormattedText(window, cur_state.stimuli(p).text, 'center', 'center',1,[],[],[],2,[],q_rects(q,:));
                        break % only one time the text has to be drawn on the screen
                    end
                    % Draw flickering faces
                    if strcmpi(cur_state.stimuli(p).type,'Texture')
                        Screen('DrawTexture', window, NewfaceTR(p),[], positions{q}(p,:),0, [], 1, [cur_freqTable(cur_state.stimuli(p).flicker,(((j-1)*12)+q)), cur_freqTable(cur_state.stimuli(p).flicker,(((j-1)*12)+q+4)),  cur_freqTable(cur_state.stimuli(p).flicker,(((j-1)*12)+q+8))]);
                    end
                end
                
                if cfg.photoDiode && strcmpi(cur_state.stimuli(1).type,'Texture') %photodiode on and face stimulus on screen
                    Screen('DrawTexture', window, diode_tex,[],diode_pos{q},0, [], 1, [cfg.diodeTable((((j-1)*12)+q)), cfg.diodeTable((((j-1)*12)+q+4)),  cfg.diodeTable((((j-1)*12)+q+8))]);
                end
                
                %eyelink fixation feedback
                if cfg.el.eyelink && ~fix && fix_loss
                    if arr_on && cfg.el.feedback
                        Screen('DrawTexture', window, arrow, [], FixArrPos{q});
                    end
                    if j-start_frame>24
                        arr_on=abs(arr_on-1); %turn on or off every 24 frames
                        start_frame=j;
                    end
                end
            end
            
        else % normal mode
            %loop over amount of things to be drawn
            for p = 1:length(cur_state.stimuli)
                if strcmpi(cur_state.stimuli(p).type,'Texture')
                    Screen('DrawTexture', window, NewfaceTR(p),[], positions(p,:),0, [], 1, [cur_freqTable(cur_state.stimuli(p).flicker,j), cur_freqTable(cur_state.stimuli(p).flicker,j),  cur_freqTable(cur_state.stimuli(p).flicker,j)]);
                end
                if strcmpi(cur_state.stimuli(p).type,'DrawLines')
                    Screen('DrawLines', window, cur_state.stimuli(p).coords{ConMat(i,1)}, cur_state.stimuli(p).params{1}, 1, cur_state.stimuli(p).params{2}, 2);
                end
                if strcmpi(cur_state.stimuli(p).type,'Text')
                    DrawFormattedText(window, cur_state.stimuli(p).text, 'center', center_nm(2),1);
                    break % only one time the text has to be drawn on the screen
                end
            end
            
            %photodiode (might not be so useful in normal mode)
            if cfg.photoDiode && strcmpi(cur_state.stimuli(1).type,'Texture') %photodiode on and face stimulus on screen
                Screen('DrawTexture', window, diode_tex,[],diode_pos,0, [], 1, [cfg.diodeTable(j), cfg.diodeTable(j),  cfg.diodeTable(j)]);
            end
            
            %Eyelink fixation feedback
            if cfg.el.eyelink && ~fix && fix_loss
                if arr_on && cfg.el.feedback
                    Screen('DrawTexture', window, arrow, [], FixArrPos_nm);
                end
                if j-start_frame>24                    
                    arr_on=abs(arr_on-1); %turn on or off every 24 frames
                    start_frame=j;
                end
            end
        end
        
        % Flip all of that to the screen
        [vbl, stim_time] = Screen('Flip', window, vbl + 0.5 * ifi);
        
         if i==1 && j==1 %first trial
            t0 = GetSecs;
            cfg.t0=t0;
            t_trial = t0;      %trial start time
            cfg.t_trial=t0;
        elseif  j==1
            t_trial = GetSecs;      %trial start time
            cfg.t_trial=t_trial;
        end
        
        %SEND TRIGGERS on state change
        if state_change
            if strcmp(cur_state.name, 'stim_onset') || strcmp(cur_state.name, 'baseline')
                cfg=sendTrigger(cfg,state_trig);
            end
            if strcmp(cur_state.name, 'movement')
                figureTime2=GetSecs; %movement onset timestamp
                figureTime3=stim_time; %movement onset timestamp derived from computed visual onset
                cfg=sendTrigger(cfg,TriggerMov);
                KbQueueFlush;%Response time starts now so clean the slate
            end
            if strcmp(cur_state.name, 'response')
                resp_start=GetSecs;
            end
            state_change=0;
        end
        
        %we want to reset the trigger 50ms after the last trigger
        if cfg.triggerSent && GetSecs>(cfg.triggerTime+0.05) && ~cfg.debugmode
            io64(cfg.ioObjTrig,cfg.PortAddress,0);
            cfg.triggerSent=0;
        end
        
        %CHECK RESPONSES
        % Check for q or e, quit or eyelink feedback off
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown
            if keyCode(cfg.el.eyelinkKey) && ~el_toggle %toggle eyelink feedback
                if cfg.el.feedback
                    cfg.el.feedback=0;
                    log(fid,'Disabling eyelink feedback');
                else
                    cfg.el.feedback=1;
                    log(fid,'Enabling eyelink feedback');
                end
                el_toggle=1;
                el_toggle_time=GetSecs;
            elseif keyCode(cfg.quitKey)
                log(fid, 'Experiment aborted');
                cleanup(cfg);
                return
            end
        end
        
        %give the toggle some time before being active again (to avoid
        %rapid switching.
        if el_toggle && GetSecs>(el_toggle_time+0.3)
            el_toggle=0;
        end
        
        while GetResp %turns at eye movement
            [keyIsDown, firstpress]=KbQueueCheck();
            
            keyCode=find(firstpress);
            
            if length(keyCode)>1 %two or more buttons pressed
                [~,ind]=min(firstpress(keyCode));
                keyCode = keyCode(ind); %select first response
            end
            
            t_keypress=firstpress(keyCode);
                        
            if ~cfg.PracticeMode
                if keyCode==cfg.KeyLeft
                    RespKey = 1;
                end
                if keyCode==cfg.KeyRight
                    RespKey = 2;
                end
            elseif keyIsDown %both keys and button boxes count
                if keyCode==cfg.KeyLeft || keyCode==cfg.BBLeft
                    RespKey = 1;
                end
                if keyCode==cfg.KeyRight || keyCode==cfg.BBRight
                    RespKey = 2;
                end
            end
                                        
            if keyIsDown
                if ismember(RespKey,[1,2])
                    rt = (t_keypress-t_trial)-figureTime;
                    rt2 = t_keypress-figureTime2;
                    rt3 = t_keypress-figureTime3;
                    
                    %send response trigger
                    cfg=sendTrigger(cfg,TriggerResp);
                    
                    %empty frame
                    Screen('Flip', window);
                    
                    %change flags
                    GetResp = 0;
                    breakFlag = 1;
                    
                    if RespKey == 1
                        if ConMat(i,1) == 1 % attend left
                            if ConMat(i,4) == 1
                                % attend left and move left
                                correct = 1;
                            elseif ConMat(i,4) == 2
                                correct = 0;
                            end
                        elseif ConMat(i,1) == 2 %attend right
                            if ConMat(i,4) == 1
                                correct = 1;
                            elseif ConMat(i,4) == 2
                                correct = 0;
                            end
                        end
                        
                    elseif RespKey == 2
                        if ConMat(i,1) == 1 % attend left
                            if ConMat(i,4) == 1
                                % attend left and move left
                                correct = 0;
                            elseif ConMat(i,4) == 2
                                correct = 1;
                            end
                        elseif ConMat(i,1) == 2 %attend right
                            if ConMat(i,4) == 1
                                correct = 0;
                            elseif ConMat(i,4) == 2
                                correct = 1;
                            end
                        end
                    end
                    
                    %Log response and RT
                    val={'incorrect','correct'};
                    key_dir={'left','right'};
                    log(fid,[int2str(cfg.block) ' ' int2str(cfg.cur_trial) ' ' num2str(t_keypress-t0)  ' ' num2str(t_keypress-t_trial) ' Response: ' val{correct+1} ' Key: ' key_dir{RespKey} ' RT1: ' num2str(rt*1000) 'ms RT2: ' num2str(rt2*1000) 'ms RT3: ' num2str(rt3*1000) 'ms']);
                    
                    if cfg.triggerSent && ~cfg.debugmode
                        WaitSecs(0.05)
                        %The trigger still has to be reset
                        io64(cfg.ioObjTrig,cfg.PortAddress,0);
                        cfg.triggerSent=0;
                    end
                    break
                else
                    if keyIsDown && keyCode==cfg.escapeKey
                        log(fid, 'Experiment aborted');
                        cleanup(cfg);
                        return
                    end
                end
            else
                
                if strcmp(cur_state.name, 'movement') % break out of while loop to go to next syncframe
                    break;
                else if GetSecs-resp_start<cfg.maxresp %check whether we are still within allowed response time
                        WaitSecs(0.001); %let's wait a whole millisecond and try again
                    else % time out
                        log(fid,[int2str(cfg.block) ' ' int2str(cfg.cur_trial) ' ' num2str(GetSecs-t0)  ' ' num2str(GetSecs-t_trial) ' MISSING Response']);
                        break
                    end
                end
            end
        end
        
        % break out of for loop if response is given
        if breakFlag
            % this break also makes sure no feedback screen will be shown.
            break
        end
    end
    
    %some logging about RT's and stuff!
    cfg.behavior(i,1) = i;           %trailNr
    cfg.behavior(i,2) = ConMat(i,1); %left or right attention
    cfg.behavior(i,3) = ConMat(i,2); %condition
    cfg.behavior(i,4) = ConMat(i,4); %left or right eye movement
    cfg.behavior(i,5) = RespKey;     %left or right key response
    cfg.behavior(i,6) = (figureframes-cueframes)*ifi;     %delay time in s
    cfg.behavior(i,7) = rt;          %in ms
    cfg.behavior(i,8) = rt2;          %in ms
    cfg.behavior(i,9) = rt3;          %in ms
    cfg.behavior(i,10) = correct;     %answer
    
    correct = NaN;
    RespKey = NaN;
    
    %Small break after every so many trials, give some feedback
    if ~mod(i,cfg.BreakTrials)
        KbQueueFlush();
        log(fid,[int2str(cfg.block) ' ' int2str(cfg.cur_trial) ' ' num2str(GetSecs-t0) ' Within-block break']);
        if i ~= nTrials
            cur_acc=(nansum(cfg.behavior(:,10))/length(cfg.behavior(:,10)))*100;
            block_acc=(nansum(cfg.behavior(i-cfg.BreakTrials+1:end,10))/length(cfg.behavior(i-cfg.BreakTrials+1:end,10)))*100;
            perf_message_pos={'Keep it up!','Great!','Nicey!','You''re doing very well!','Good stuff.','Awesome!','Sweet.','Groovy','Doing great','That''s it!','Splendid :)'};
            perf_message_med={'Good effort, try to stay focused.','Well done, but try to improve your score.','Not bad at all, but stay on it!','Good! You missed some though.'};
            perf_message_lo={'Keep trying.','Try to improve your score a bit','Don''t give up!','Oops, not too great, make sure to pay attention to the eye movement direction.','It''s ok to make mistakes, but keep trying.'};
            if block_acc>75
                perf_mess=perf_message_pos{round(rand*(length(perf_message_pos)-1))+1};
            else if block_acc>60
                    perf_mess=perf_message_med{round(rand*(length(perf_message_med)-1))+1};
                else
                    perf_mess=perf_message_lo{round(rand*(length(perf_message_lo)-1))+1};
                end
            end
            cur_rt=nanmean(cfg.behavior(i-cfg.BreakTrials+1:end,8));
            if cur_rt>0.8 && block_acc>60
                rt_add='However, try to respond a bit faster.';
            end
            fb_message=sprintf('Take a small Break \nYou have completed %i of %i trials in this block \nYour accuracy this block was %i%%. %s\nYour overall accuracy is %i%%.\nPress any Key to continue',i,nTrials,round(block_acc),perf_mess,round(cur_acc));
            if cfg.rapidMode
                for q = 1 : 4
                    [~,~,~]=DrawFormattedText(window,fb_message, 'center', 'center',1,[],[],[],2,[],q_rects(q,:));
                end
            else
                DrawFormattedText(window, fb_message, 'center', center_nm(2),1);
            end
            
            vbl = Screen('Flip', window);
            WaitSecs(1);
            
            GetResp=1;
            while GetResp %turns at eye movement
                [keyIsDown, ~]=KbQueueCheck();
                if keyIsDown
                    GetResp=0;
                end
            end
            
            %add empty screen
            Screen('Flip', window);
            WaitSecs(0.1);
        end
    end
end

%End of block message
if cfg.rapidMode
    for q = 1 : 4
        [~,~,~]=DrawFormattedText(window,'End of this block! \nWell done' , 'center', 'center',1,[],[],[],2,[],q_rects(q,:));
    end
else
    DrawFormattedText(window,'End of this block \nWell Done!', 'center', center_nm(2),1);
end

vbl = Screen('Flip', window);

log(fid, ['END of Block ' int2str(cfg.block)]);

%stop eyelink & transfer file
if cfg.el.eyelink
    Eyelink('Message', 'END OF SESSION');
    el_Stop(cfg);
end

%release KbQueue
KbQueueRelease();


%return to lower priority
if ~cfg.debugmode
    Priority(0);
end

%close logfile
fclose(fid);

%save cfg
save([cfg.exp_dir 'Logs' filesep cfg.sub filesep 'cfg_' cfg.sub '_ses_' int2str(cfg.block) '_' cfg.datestamp '.mat'],'cfg');

WaitSecs(2);

%set propixx to normal state
if cfg.rapidMode && ~cfg.debugmode || cfg.DataPixxOnly
    Datapixx('SetPropixxDlpSequenceProgram', 0);
    Datapixx('RegWrRd');
    Datapixx('close');
end

%Close screen
sca

disp('Experimental session finished');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SubFunctions

% function to create log files
function [] = log(fid,msg)
fprintf(fid,[msg '\n']);
disp(msg)
end

%function to send MEG and eyelink triggers, and log triggers in logfile
function [cfg]= sendTrigger(cfg,trig)

%send trigger to MEG
if ~cfg.debugmode
    io64(cfg.ioObjTrig,cfg.PortAddress,trig);
    cfg.triggerTime=GetSecs;
    cfg.triggerSent=1;
else
    cfg.triggerTime=GetSecs;
end

%send trigger to eyelink
if cfg.el.eyelink
    Eyelink('Message', [int2str(cfg.block) ' ' int2str(cfg.cur_trial) ' ' num2str(cfg.triggerTime-cfg.t0) ' ' num2str(cfg.triggerTime-cfg.t_trial) ' Trigger: ' int2str(trig)]);
end

%Log trigger
fprintf(cfg.fid,[int2str(cfg.block) ' ' int2str(cfg.cur_trial) ' ' num2str(cfg.triggerTime-cfg.t0) ' ' num2str(cfg.triggerTime-cfg.t_trial)  ' Trigger: ' int2str(trig) '\n']);

%log triggers in cfg
cfg.Triggers=[cfg.Triggers ; cfg.block cfg.cur_trial cfg.triggerTime-cfg.t0 cfg.triggerTime-cfg.t_trial trig];

end

function [] = cleanup(cfg)

log(cfg.fid, ['END of Block ' int2str(cfg.block) ' -ABORTED-']);

%release KbQueue
KbQueueRelease();


%Return propixx to normal state
if cfg.rapidMode && ~cfg.debugmode
    Datapixx('SetPropixxDlpSequenceProgram', 0);
    Datapixx('RegWrRd');
    Datapixx('close');
end

%lower priority
if ~cfg.debugmode
    Priority(0);
end


%stop eyelink
if cfg.el.eyelink
    Eyelink('Message', 'END OF SESSION - ABORTED');
    el_Stop(cfg);
end

%close logfile
fclose(cfg.fid);

%save cfg
save([cfg.exp_dir 'Logs' filesep cfg.sub filesep 'cfg_' cfg.sub '_ses_' int2str(cfg.block) '_' cfg.datestamp '.mat'],'cfg');

%close screen
sca

%throw warning due to prematurely aborted experiment
warning('Experiment aborted');
end

