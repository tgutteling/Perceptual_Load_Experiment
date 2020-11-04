function []=run_intro(cfg)
% Intro
message=cell(1,10);
drawing=cell(1,10);
texture=cell(1,10);
message_yshift(1:10)=85;%defaults shift

message{1}='Welcome to the experiment!\nThroughout the experiment, please keep your eyes\non the fixation cross.\nAlso, make sure you''re not accidentally pressing any buttons\nPress any button to continue...';
drawing{1}={'DrawLines', cfg.window, cfg.CrossCoords, cfg.lineWidthPix, 1};

message{2}='Every trial you will see two faces,\nOne on the left and one on the right\n(Press any button to continue...)';

texture{2}{1}={'DrawTexture', cfg.window,cfg.faceTex(6,4),[]};
texture{2}{2}={'DrawTexture', cfg.window,cfg.faceTex(6,4),[]};

message{3}='An arrow will indicate which face you have to pay attention to\n but you have to keep looking at the fixation cross!\n(Press any button to continue...)';

drawing{3}={'DrawLines', cfg.window, cfg.LeftCueCoords, cfg.lineWidthPix, 1};
texture{3}{1}={'DrawTexture', cfg.window,cfg.faceTex(6,4),[]};
texture{3}{2}={'DrawTexture', cfg.window,cfg.faceTex(6,4),[]};

message{4}='After some time, both faces will make an eye movement\nThe faces will disappear shorty after this\nYour task is to indicate whether the selected\n face (left in this case) looked left or right\n(Press any button to continue...)';
message_yshift(4)=95;
drawing{4}={'DrawLines', cfg.window, cfg.CrossCoords, cfg.lineWidthPix, 1};
texture{4}{1}={'DrawTexture', cfg.window,cfg.faceTex(6,4),[]};
texture{4}{2}={'DrawTexture', cfg.window,cfg.faceTex(6,4),[]};

message{5}='Press the left button (index finger) when you see\nthe eyes moving to the left\n(Press any button to continue...)';
drawing{5}={'DrawLines', cfg.window, cfg.CrossCoords, cfg.lineWidthPix, 1};
texture{5}{1}={'DrawTexture', cfg.window,cfg.faceTex(6,5),[]};
texture{5}{2}={'DrawTexture', cfg.window,cfg.faceTex(6,5),[]};

message{6}='..and press the right button when you see\nthe eyes moving to the right\n(Press any button to continue...)';
drawing{6}={'DrawLines', cfg.window, cfg.CrossCoords, cfg.lineWidthPix, 1};
texture{6}{1}={'DrawTexture', cfg.window,cfg.faceTex(6,6),[]};
texture{6}{2}={'DrawTexture', cfg.window,cfg.faceTex(6,6),[]};

message{7}='Sometimes the left and right face look in different directions\nRemember to pay attention to the selected face\nTry to be both fast and accurate\nand stay within the time limit (1 second)\n(Press any button to continue...)';
drawing{7}={'DrawLines', cfg.window, cfg.CrossCoords, cfg.lineWidthPix, 1};
texture{7}{1}={'DrawTexture', cfg.window,cfg.faceTex(6,6),[]};
texture{7}{2}={'DrawTexture', cfg.window,cfg.faceTex(6,5),[]};

message{8}='Sometimes one or both faces will be noisy, this\nmakes it a bit harder, but try your best!\nYou might see the faces flickering a little bit, this is normal\n(Press any button to continue...)';
drawing{8}={'DrawLines', cfg.window, cfg.CrossCoords, cfg.lineWidthPix, 1};
texture{8}{1}={'DrawTexture', cfg.window,cfg.faceTex(6,2),[]};
texture{8}{2}={'DrawTexture', cfg.window,cfg.faceTex(6,5),[]};

message{9}='Try to keep your eyes on the fixation cross\n and try not to blink too much\n(you have to blink at some point of course)\nIt is best to blink in between trials, when the fixation cross disappears\n(Press any button to continue...)';

message{10}='There will be two blocks of 20 minutes with a break in between,\nand smaller breaks every 2-3 minutes\nIf anything is unclear, let the experimenter know\nOtherwise, GOOD LUCK!\n(Press any button to continue...)';

%set up KbQueue
KbQueueCreate;
KbQueueStart(); %Start listening
KbQueueFlush();

for m=1:length(message)
    
    KbQueueFlush();
    %Messages
    if cfg.rapidMode
        for q=1:4
            [~,~,~]=DrawFormattedText(cfg.window, message{m}, 'center', 'center',1,[],[],[],[],[],cfg.q_rects(q,:)-[0 0 0 message_yshift(m)]);
            if ~isempty(drawing{m})
                Screen(drawing{m}{:},cfg.qcenters(q,:),2)
            end
            if ~isempty(texture{m})
                for p=1:2
                    Screen(texture{m}{p}{:},cfg.positions{q}(p,:),0, [], 1, [0.5 0.5 0.5]);
                end
            end
        end
    else
        DrawFormattedText(cfg.window, message{m}, 'center', cfg.center_nm(2)-message_yshift(m)*2,1)
        if ~isempty(drawing{m})
            Screen(drawing{m}{:},cfg.center_nm,2)
        end
        if ~isempty(texture{m})
            for p=1:2
                Screen(texture{m}{p}{:},cfg.positions(p,:),0, [], 1, [0.5 0.5 0.5]);
            end
        end
    end
    vbl = Screen('Flip', cfg.window);
    
    WaitSecs(.4);
    
    GetResp=1;
    while GetResp 
        [keyIsDown, ~]=KbQueueCheck();
        if keyIsDown
            GetResp=0;
        end
    end
    %KbQueueWait;
end

KbQueueRelease();

end

