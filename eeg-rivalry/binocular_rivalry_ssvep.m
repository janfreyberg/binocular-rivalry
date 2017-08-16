
function binocular_ssvep

clearvars;

global pxsize frameWidth ycen xcen fixWidth scr l_key u_key d_key r_key...
    esc_key stimRect fixLines fixPoint frameRect checkerboard_arrays...
    checkerboard_alpha frequencies period_frames trialdur framedur...
    breakdur currTrial

try
% get basic info, set filename.
rng('shuffle');
% subject info and screen info

subject_info = inputdlg({'Subject ID', 'Today''s date'}, 'Subject Info', 1, {'test', date});
scr_diagonal = 24;
scr_distance = 60;

savefile = fullfile(pwd, 'data', [date, '-', subject_info{1}, '.mat']);

% experiment Variables.
scr_background = 127.5;
scr_no = 0;
scr_dimensions = Screen('Rect', scr_no);
xcen = scr_dimensions(3)/2;
ycen = scr_dimensions(4)/2;

% Frame Duration
framedur = 1/120;

% Frequencies
frequencies = [5.0, 15];
period_seconds = 1 ./ frequencies;
period_frames = period_seconds / framedur;

% Trialtime in seconds
trialdur = 5;
breakdur = 2;

% stimsize in degree
stimsize = 6;

% set up Keyboard, Screen, Sound
% Keyboard
KbName('UnifyKeyNames');
u_key = KbName('UpArrow');
d_key = KbName('DownArrow');
l_key = KbName('LeftArrow');
r_key = KbName('RightArrow');
esc_key = KbName('Escape');
space_key = KbName('Space');
ent_key = KbName('Return'); ent_key = ent_key(1);
keyList = zeros(1, 256);
keyList([u_key, d_key, l_key, r_key, esc_key, ent_key, space_key]) = 1;
KbQueueCreate([], keyList); clear keyList
ListenChar(2);

% I/O driver
% config_io;
address = hex2dec('DFB8');
% in the triggers, 0=no, 1=right, 10=up, 11=up+right, 100=left,
% 101=left+right, 110=left+up, 111=left+up+right

% Open Window
Screen('Preference', 'SkipSyncTests', 1);
scr = Screen('OpenWindow', scr_no, scr_background);
HideCursor;
Screen('BlendFunction', scr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% prepare stimuli
% Stimsize
pxsize = 500;%visual_angle2pixel(6, scr_diagonal, scr_distance, scr_no);

% Make Stimuli
% This will be a 2 x 2 cell
% first row is more spokes
[checkerboard_arrays{1, 1}, checkerboard_arrays{1, 2}] = make_circular_checkerboard_pattern(4, 8, pxsize);
% second row is more circles
[checkerboard_arrays{2, 1}, checkerboard_arrays{2, 2}] = make_circular_checkerboard_pattern(8, 4, pxsize);
% alpha will be a circle
checkerboard_alpha = 255 * Circle(pxsize/2);


% Vergence Cues
fixWidth = visual_angle2pixel(stimsize / 40, scr_diagonal, scr_distance, scr_no);
fixLength = visual_angle2pixel(stimsize / 15, scr_diagonal, scr_distance, scr_no);
frameWidth = visual_angle2pixel( stimsize / 30, scr_diagonal, scr_distance, scr_no);
fixLines = [-fixLength, +fixLength, 0, 0; 0, 0, -fixLength, +fixLength];



% Demonstrate

demonstrate;

trialbreak;

% make schedule
j = 0;
for iflicker = 0:2
   for icolor = 0:1
        for itrialNo = 1:2
            j = j+1;
            if iflicker==2
                frequency_order(1:2, j) = 0;
            else
                frequency_order(1:2, j) = [iflicker, ~iflicker]+1;
            end
            color_order(1:2, j) = [icolor, ~icolor]+1;
        end
   end
end
sched = randperm(j);
frequency_order(1:2, :) = frequency_order(1:2, sched);
color_order(1:2, :) = color_order(1:2, sched);


% run trials
timestamps = cell(j, 1);
responses = cell(j, 1);
for currTrial = 1:j
    % run the trial
    [timestamps{currTrial}, responses{currTrial}] = ...
        trial(color_order(1:2, currTrial), frequency_order(1:2, currTrial));
    % have a break
    trialbreak;
end



    % Clean & close
    KbQueueFlush;
    KbQueueStop;
    KbQueueRelease;
    ListenChar(0);
    sca;
%     save(savefile);
    save('temp_binocular_ssvep.mat');

catch err
    KbQueueFlush;
    KbQueueStop;
    KbQueueRelease;
    ListenChar(0);
    sca;
%     save(savefile);
    save('temp_binocular_ssvep.mat');
    rethrow(err);
end
end


%% FUNCTIONS
function pressed = checkkey()
    global esc_key
    [~, pressed] = KbQueueCheck;
    if pressed(esc_key)
        error('Interrupted in the break!');
    end
end

function pressed = waitkey()
    global esc_key
    [~, pressed] = KbStrokeWait;
    if pressed(esc_key)
        error('Interrupted in the break!');
    end
end

function trigger(value)
    global address
%     outp(address, 0);
%     outp(address, round(value));
    WaitSecs(0.002);
%     outp(address, 0);
    disp(value);
end

function [timestamps, buttons] = trial(color_order, frequency_order)
    global pxsize frameWidth ycen xcen fixWidth scr l_key u_key d_key r_key...
    esc_key stimRect fixLines fixPoint frameRect checkerboard_arrays...
    checkerboard_alpha frequencies period_frames trialdur framedur currTrial

    % Make the textures
    % stim
    for stim1= 1:2
        for stim2= 1:2
            tmparray = cat(3, 255*checkerboard_arrays{color_order(1), stim1},...
                                zeros(size(checkerboard_arrays{stim1})),...
                                255*checkerboard_arrays{color_order(2), stim2},...
                                checkerboard_alpha);
            textures{stim1, stim2} = Screen('MakeTexture', scr, tmparray);
        end
    end
    % fixation
    fixpatch = cat(3, repmat(127.5*Circle(round(pxsize ./ 8)), 1, 1, 3), 255*Circle(round(pxsize ./ 8)));
    fixpatch = Screen('MakeTexture', scr, fixpatch);
    
    
    % start trial
    DrawFormattedText(scr, ['Trial ' num2str(currTrial) ' Ready?',...
                            '\n\nPress and hold [UP] to start.'], 'center', ycen-pxsize/3, 0);
    Screen('DrawLines', scr, fixLines, fixWidth, 0, [xcen, ycen]);
    Screen('Flip', scr);
    waitkey;
    
    
    
    % how many frames?
    trialframes = trialdur / framedur;
    timestamps = zeros(trialframes, 1);
    buttons = zeros(trialframes, 3);
    
    Priority(1);
    KbQueueStart;
    t0 = GetSecs;
    for i = 2:trialframes
        if all(frequency_order)
            stim1 = square(i*2*pi/period_frames(frequency_order(1)))/2+1.5;
            stim2 = square(i*2*pi/period_frames(frequency_order(2)))/2+1.5;
        else
            stim1 = 1;
            stim2 = 1;
        end
        
        Screen('DrawTexture', scr, textures{stim1, stim2});
        Screen('DrawTexture', scr, fixpatch);
        Screen('DrawLines', scr, fixLines, fixWidth, 0, [xcen, ycen]);
        timestamps(i) = Screen('Flip', scr) - t0;
        resp = checkkey;
        buttons(i, :) = logical(resp([l_key, u_key, r_key]));
        if any(buttons(i, :) ~= buttons(i-1, :)) && any(resp)
            trigger(bin2dec([num2str(buttons(i, 1)), num2str(buttons(i, 2)), num2str(buttons(i, 3))]) + 1);
        end
    end
    
    % clean up
    KbQueueStop;
    Priority(0);
    Screen('Close', [textures{:}]);
    Screen('DrawLines', scr, fixLines, fixWidth, 0, [xcen, ycen]);
    Screen('Flip', scr);
end

function demonstrate()
    global pxsize frameWidth ycen xcen fixWidth scr l_key u_key d_key r_key...
    esc_key stimRect fixLines fixPoint frameRect checkerboard_arrays...
    checkerboard_alpha frequencies period_frames
    

    fixpatch = cat(3, repmat(127.5*Circle(round(pxsize ./ 8)), 1, 1, 3), 255*Circle(round(pxsize ./ 8)));
    fixpatch = Screen('MakeTexture', scr, fixpatch);
    
    % Make the red grating
    colors = [1, 2];
    for stim1= 1:2
        for stim2= 1:2
            tmparray = cat(3, 255*checkerboard_arrays{colors(1), stim1},...
                                zeros(size(checkerboard_arrays{stim1})),...
                                zeros(size(checkerboard_arrays{stim1})),...
                                checkerboard_alpha);
            textures{stim1, stim2} = Screen('MakeTexture', scr, tmparray);
        end
    end
    % Flicker the red
    i = 1;
    KbQueueStart;
    while true
        Screen('DrawTexture', scr, textures{square(i*2*pi/period_frames(1))/2+1.5, square(i*2*pi/period_frames(2))/2+1.5});
        Screen('DrawTexture', scr, fixpatch);
        Screen('DrawLines', scr, fixLines, fixWidth, 0, [xcen, ycen]);
        Screen('Flip', scr);
        if any(checkkey)
            break
        end
        i = i+1;
    end
    
    % Make the blue grating
    colors = [1, 2];
    for stim1= 1:2
        for stim2= 1:2
            tmparray = cat(3, zeros(size(checkerboard_arrays{stim1})),...
                                zeros(size(checkerboard_arrays{stim1})),...
                                255*checkerboard_arrays{colors(2), stim2},...
                                checkerboard_alpha);
            textures{stim1, stim2} = Screen('MakeTexture', scr, tmparray);
        end
    end
    % Flicker the blue
    i = 1;
    KbQueueStart;
    while true
        Screen('DrawTexture', scr, textures{square(i*2*pi/period_frames(1))/2+1.5, square(i*2*pi/period_frames(2))/2+1.5});
        Screen('DrawTexture', scr, fixpatch);
        Screen('DrawLines', scr, fixLines, fixWidth, 0, [xcen, ycen]);
        Screen('Flip', scr);
        if any(checkkey)
            break
        end
        i = i+1;
    end
    
    % Make the blue grating
    colors = [1, 2];
    for stim1= 1:2
        for stim2= 1:2
            tmparray = cat(3, 255*checkerboard_arrays{colors(1), stim1},...
                                zeros(size(checkerboard_arrays{stim1})),...
                                255*checkerboard_arrays{colors(2), stim2},...
                                checkerboard_alpha);
            textures{stim1, stim2} = Screen('MakeTexture', scr, tmparray);
        end
    end
    % Flicker the blue
    i = 1;
    KbQueueStart;
    while true
        Screen('DrawTexture', scr, textures{square(i*2*pi/period_frames(1))/2+1.5, square(i*2*pi/period_frames(2))/2+1.5});
        Screen('DrawTexture', scr, fixpatch);
        Screen('DrawLines', scr, fixLines, fixWidth, 0, [xcen, ycen]);
        Screen('Flip', scr);
        if any(checkkey)
            break
        end
        i = i+1;
    end
    
    Screen('Close', [textures{:}]);
end



function trialbreak()
    % Take a Break
    global scr breakdur ycen pxsize
    KbQueueStart;
    for lapsedTime = 0:breakdur
        DrawFormattedText(scr, ['Break for ' num2str(breakdur-lapsedTime)], 'center', ycen-pxsize/3, 0);
        Screen('Flip', scr);
        checkkey;
        WaitSecs(1);
    end
    KbQueueStop;
end

