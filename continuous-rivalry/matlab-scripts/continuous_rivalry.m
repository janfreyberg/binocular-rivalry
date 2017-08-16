function rivalryblocks
try
    Screen('Preference', 'SkipSyncTests', 1);
%% Session Variables

ID = input('Participant ID? ', 's');
scr_diagonal = input('Screen Diagonal? ');
scr_distance = 57;
diagnosis = input('Diagnosis? ');

tstamp = clock;
if ~isdir( fullfile(pwd, 'Results', num2str(diagnosis)) )
    mkdir( fullfile(pwd, 'Results', num2str(diagnosis)) );
end
savefile = fullfile(pwd, 'Results', num2str(diagnosis), [sprintf('-%02d-%02d-%02d-%02d%02d-', tstamp(1), tstamp(2), tstamp(3), tstamp(4), tstamp(5)), ID, '.mat']);



%% Experiment Variables

scr_background = 127.5;
scr_no = max(Screen('Screens'));
scr_dimensions = Screen('Rect', scr_no);
xcen = scr_dimensions(3)/2;
ycen = scr_dimensions(4)/2;

freq = 6;

stimsizedeg = 3.5;
stimsize = visual_angle2pixel(stimsizedeg, scr_diagonal, scr_distance, scr_no);

gratingsizedeg = 3.5;
gratingsize = visual_angle2pixel(gratingsizedeg, scr_diagonal, scr_distance, scr_no);

blockdur = 40;
samplerate = 1/(60*4); % time in between sample collection

iti = 2;

reps = 6; % repitions per condition



%% Set-Up

% Keyboard
KbName('UnifyKeyNames');
l_key = KbName('LeftArrow'); r_key = KbName('RightArrow');
u_key = KbName('UpArrow'); d_key = KbName('DownArrow');
esc_key = KbName('Escape');
ent_key = KbName('Return'); ent_key = ent_key(1);
ListenChar(2);

% Sound
InitializePsychSound;
pa = PsychPortAudio('Open', [], [], [], [], [], 256);
bp400 = PsychPortAudio('CreateBuffer', pa, [MakeBeep(400, 0.2); MakeBeep(400, 0.2)]);
PsychPortAudio('FillBuffer', pa, bp400);

% Screen
scr = Screen('OpenWindow', scr_no, scr_background);
HideCursor;
if IsLinux
    Screen('TextFont', scr, '-schumacher-clean-medium-r-normal--0-0-75-75-c-0-koi8-r');
end
Screen('BlendFunction', scr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Priority
Priority(1);



%% Stimuli

% Frame
frameWidth = 5;
frameColour = 0;

% Fixation
fixWidth = 4;
fixLength = 8;
fixColour = 0;

% Grating
freq = 3;
contrast = 0.8;



%% Find the spot

% this fct only presents two lines. try and make them match up with mirrors
calibrate_mirrors;

[offset, frameRect, stimRect, fixLines, gratRect] = find_offset;

WaitSecs(1);
Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
Screen('Flip', scr);
WaitSecs(0.5);
KbWait;



%% Demonstrate Stimuli
% textures(1) = make_image(fullfile(pwd, 'Photos', '7.jpg'), [255 0 0]);
% textures(2) = make_image(fullfile(pwd, 'Photos', '10.jpg'), [0 255 0]);
% 
% demonstrate_stimuli(textures, 'images');
% Screen('Close', textures);
% 
% textures(1) = make_grating(45, 1);
% textures(2) = make_grating(-45, 1);
% 
% demonstrate_stimuli(textures, 'gratings');
% Screen('Close', textures);
% 
% textures(1) = make_image(fullfile(pwd, 'Photos', '1.jpg'), [255 0 0]);
% textures(2) = make_image(fullfile(pwd, 'Photos', '2.jpg'), [0 255 0]);
% 
% rivalry_block(textures, 15, 'images');
% Screen('Close', textures);
% take_break;
% 
% textures(1) = make_grating(45, 1);
% textures(2) = make_grating(-45, 1);
% 
% practice = rivalry_block(textures, 15, 'gratings');
% Screen('Close', textures);
% take_break;


%% Blocks of Rivalry

hj = 0;
for h = 1:2
   for j = 1:3
       for iLR = 1:2
           for t = 1:reps/2
            hj = hj + 1;
            stimType(hj) = h;
            trialType(hj) = j;
            LR(hj) = iLR;
           end
       end
   end
end

order = randperm(hj);
stimType = stimType(order);
trialType = trialType(order);
LR = LR(order);

trialType = [3 3];
stimType = [1 2];

for h = 1:hj
    
    j = find(find(stimType==stimType(h) & trialType == trialType(h))==h);
    
    if trialType(h) == 1 && stimType(h) == 1
        %% Blocked Rivalry Colour
        block(trialType(h), stimType(h)).type = 'images';

        block(trialType(h), stimType(h)).imgNrs(j, 1:2) = randsample(15, 2)';

        image1 = fullfile(pwd, 'Photos', [num2str(block(trialType(h), stimType(h)).imgNrs(j, 1)), '.jpg']);
        image2 = fullfile(pwd, 'Photos', [num2str(block(trialType(h), stimType(h)).imgNrs(j, 2)), '.jpg']);

        textures(LR(h)) = make_image(image1, [255, 0, 0]);
        textures(mod(LR(h), 2)+1) = make_image(image2, [0, 255, 0]);

        [block(trialType(h), stimType(h)).pressSecs(:,:,j), block(trialType(h), stimType(h)).pressList(:,:,j)] = rivalry_block(textures, blockdur, block(trialType(h), stimType(h)).type);

        Screen('Close', textures);
        if h < hj
            take_break;
        elseif h == round(hj/2)
            take_long_break;
        end
    
    elseif trialType(h) == 1 && stimType(h) == 2
        %% Blocked Rivalry Gratings
        block(trialType(h), stimType(h)).type = 'gratings';

        textures(LR(h)) = make_grating( 45, 1 );
        textures(mod(LR(h), 2)+1) = make_grating( -45, 1 );

        [block(trialType(h), stimType(h)).pressSecs(:,:,j), block(trialType(h), stimType(h)).pressList(:,:,j)] = rivalry_block(textures, blockdur, block(trialType(h), stimType(h)).type);

        Screen('Close', textures);
        if h < hj
            take_break;
        elseif h == round(hj/2)
            take_long_break;
        end
    
    elseif trialType(h) == 2 && stimType(h) == 1
        %% Simulation - Images
        block(trialType(h), stimType(h)).type = 'images';
        block(trialType(h), stimType(h)).imgNrs(j, 1:2) = randsample(15, 2)';
        

        [block(trialType(h), stimType(h)).pressSecs(:,:,j), block(trialType(h), stimType(h)).pressList(:,:,j), block(trialType(h), stimType(h)).course(:,:,j)] = rivalry_simulation(block(trialType(h), stimType(h)).type, blockdur, LR(h), block(trialType(h), stimType(h)).imgNrs(j, 1:2));
        if h < hj
            take_break;
        elseif h == round(hj/2)
            take_long_break;
        end
        
    elseif trialType(h) == 2 && stimType(h) == 2
        %% Simulation - Gratings
        block(trialType(h), stimType(h)).type = 'gratings';
        

        [block(trialType(h), stimType(h)).pressSecs(:,:,j), block(trialType(h), stimType(h)).pressList(:,:,j), block(trialType(h), stimType(h)).course(:,:,j)] = rivalry_simulation(block(trialType(h), stimType(h)).type, blockdur, LR(h));
        if h < hj
            take_break;
        elseif h == round(hj/2)
            take_long_break;
        end

    elseif trialType(h) == 3 && stimType(h) == 1
        %% Simulation - Images
        block(trialType(h), stimType(h)).type = 'images';
        block(trialType(h), stimType(h)).imgNrs(j, 1:2) = randsample(15, 2)';
        

        [block(trialType(h), stimType(h)).pressSecs(:,:,j), block(trialType(h), stimType(h)).pressList(:,:,j), block(trialType(h), stimType(h)).course(:,:,j)] = rivalry_simulation(block(trialType(h), stimType(h)).type, blockdur, LR(h), block(trialType(h), stimType(h)).imgNrs(j, 1:2), 'sudden');
        if h < hj
            take_break;
        elseif h == round(hj/2)
            take_long_break;
        end
        
    elseif trialType(h) == 3 && stimType(h) == 2
        %% Simulation - Gratings
        block(trialType(h), stimType(h)).type = 'gratings';
        

        [block(trialType(h), stimType(h)).pressSecs(:,:,j), block(trialType(h), stimType(h)).pressList(:,:,j), block(trialType(h), stimType(h)).course(:,:,j)] = rivalry_simulation(block(trialType(h), stimType(h)).type, blockdur, LR(h), [], 'sudden');
        if h < hj
            take_break;
        elseif h == round(hj/2)
            take_long_break;
        end
        
    end
end



%% Shutdown
sca;
save(savefile);
ListenChar(0);
if strcmp(input('Do you want to keep the data? y / n >> ', 's'), 'n')
    delete(savefile);
    disp('Data not saved.');
end
Priority(0);



catch err
%% Catch
sca;
savefile = [savefile(1:(size(savefile, 2)-4)), '-ERROR.mat'];
save(savefile);
ListenChar(0);
if strcmp(input('Do you want to keep the data? y / n >> ', 's'), 'n')
    delete(savefile);
    disp('Data not saved.');
end
Priority(0);
rethrow(err);



end
%% Calibration, Break and Stimulus Creation

    function [offset, frameRect, stimRect, fixLines, gratRect] = find_offset
        offset = stimsize;
        % Firstly, present white & black circles to avoid merging and find
        % the spot where they meet!
        while 1
            frameRect= [-sqrt(0.5 * stimsize^2)-frameWidth;
                        -sqrt(0.5 * stimsize^2)-frameWidth;
                         sqrt(0.5 * stimsize^2)+frameWidth;
                         sqrt(0.5 * stimsize^2)+frameWidth];
            frameRect= [xcen-offset, xcen+offset; ycen, ycen; xcen-offset, xcen+offset; ycen, ycen] + [frameRect, frameRect];
            stimRect= [-stimsize/2;
                       -stimsize/2;
                        stimsize/2;
                        stimsize/2];
            stimRect= [xcen-offset, xcen+offset; ycen, ycen; xcen-offset, xcen+offset; ycen, ycen] + [stimRect, stimRect];
            fixLines = [-fixLength, +fixLength, 0, 0;
                        0, 0, -fixLength, +fixLength];
            fixLines = [xcen-offset, xcen-offset, xcen-offset, xcen-offset, xcen+offset, xcen+offset, xcen+offset, xcen+offset; ycen, ycen, ycen, ycen, ycen, ycen, ycen, ycen]+...
                        [fixLines, fixLines];
            
                        
            Screen('FrameRect', scr, [255,0;255,0;255,0], frameRect, frameWidth);
            Screen('DrawLines', scr, fixLines, fixWidth, [255,255,255,255,0,0,0,0;255,255,255,255,0,0,0,0;255,255,255,255,0,0,0,0]);
            Screen('Flip', scr);
            
            WaitSecs(0.2);
            
            [~, keyCode] = KbWait;
            
            if keyCode(u_key)
                offset = offset + 10;
            elseif keyCode(d_key)
                offset = offset - 10;
            elseif keyCode(l_key)
                offset = offset + 2;
            elseif keyCode(r_key)
                offset = offset - 2;
            elseif keyCode(ent_key)
                offset = offset + stimsize/2 + frameWidth;
                break;
            elseif keyCode(esc_key)
                error('You interrupted the script!');
            end
            if offset < stimsize
                offset = stimsize;
            end
        end
        % Once the point where they meet has been established, make a jump
        % and then present the circles the same colour. Also present text
        % to verify
        while 1
            frameRect= [-sqrt(0.5 * stimsize^2)-frameWidth;
                        -sqrt(0.5 * stimsize^2)-frameWidth;
                         sqrt(0.5 * stimsize^2)+frameWidth;
                         sqrt(0.5 * stimsize^2)+frameWidth];
            frameRect= [xcen-offset, xcen+offset; ycen, ycen; xcen-offset, xcen+offset; ycen, ycen] + [frameRect, frameRect];
            stimRect= [-stimsize/2;
                       -stimsize/2;
                        stimsize/2;
                        stimsize/2;];
            stimRect= [xcen-offset, xcen+offset; ycen, ycen; xcen-offset, xcen+offset; ycen, ycen] + [stimRect, stimRect];
            fixLines = [-fixLength, +fixLength, 0, 0;
                        0, 0, -fixLength, +fixLength];
            fixLines = [xcen-offset, xcen-offset, xcen-offset, xcen-offset, xcen+offset, xcen+offset, xcen+offset, xcen+offset; ycen, ycen, ycen, ycen, ycen, ycen, ycen, ycen]+...
                        [fixLines, fixLines];
            
            
            [newX, ~] = Screen('DrawText', scr, 'A', xcen-offset-stimsize/3, ycen-stimsize/3, [255 0 0]);
            width = newX - (xcen-offset-stimsize/3); height = Screen('TextSize', scr);
            Screen('DrawText', scr, 'B', xcen-offset+stimsize/3-width, ycen-stimsize/3, [0 255 0]);
            Screen('DrawText', scr, 'C', xcen-offset-stimsize/3, ycen+stimsize/3-height-2, [0 0 255]);
            Screen('DrawText', scr, 'D', xcen-offset+stimsize/3-width, ycen+stimsize/3-height-2, [0 255 255]);
            
            Screen('DrawText', scr, 'A', xcen+offset-stimsize/3, ycen-stimsize/3, [255 0 0]);
            Screen('DrawText', scr, 'B', xcen+offset+stimsize/3-width, ycen-stimsize/3, [0 255 0]);
            Screen('DrawText', scr, 'C', xcen+offset-stimsize/3, ycen+stimsize/3-height-2, [0 0 255]);
            Screen('DrawText', scr, 'D', xcen+offset+stimsize/3-width, ycen+stimsize/3-height-2, [0 255 255]);
            
            Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
            Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
            Screen('Flip', scr);
            
            WaitSecs(0.2);
            
            [~, keyCode] = KbWait;
            
            if keyCode(u_key)
                offset = offset + 10;
            elseif keyCode(d_key)
                offset = offset - 10;
            elseif keyCode(l_key)
                offset = offset + 2;
            elseif keyCode(r_key)
                offset = offset - 2;
            elseif keyCode(ent_key)
                break;
            elseif keyCode(esc_key)
                error('You interrupted the script!');
            end
        end
        gratRect = stimRect + [(stimsize-gratingsize)/2, (stimsize-gratingsize)/2;
                        (stimsize-gratingsize)/2, (stimsize-gratingsize)/2;
                        -(stimsize-gratingsize)/2, -(stimsize-gratingsize)/2;
                        -(stimsize-gratingsize)/2, -(stimsize-gratingsize)/2];
    end


    function calibrate_mirrors
        % present two lines, one blue, one red (prevent fusing!), and make
        % sure the vertical alignment is calibrated properly
        Screen('DrawLine', scr, [255, 0, 0], 0, ycen, xcen, ycen, 5);
        Screen('DrawLine', scr, [0, 0, 255], xcen, ycen, 2*xcen, ycen, 5);
        Screen('Flip', scr);
        WaitSecs(2);
        KbWait;
        WaitSecs(2);
    end



    function tex = make_grating(angle, contrast, alpha)
        
        [x, y] = meshgrid(linspace(-1.6/2, 1.6/2, stimsize));
        
        
        
        rings = cos(2*pi*sqrt(x.^2 + y.^2)/1.2 - pi/1.2);
        rings( sqrt(x.^2 + y.^2) < 0.5 ) = 1;
        rings( rings<0 ) = 0;
        
        
        x2 = x * cosd(angle); y2 = y * sind(angle);
        wave = 127.5+ contrast* 127.5* cos((x2 + y2)* 2*pi * 3.25 ) .*rings;
        
        if nargin < 3 || isempty(alpha)
            alpha = 255*ones(size(wave));
        end
        
        
        wave = cat(3, wave, wave, wave, alpha);
        
        tex = Screen('MakeTexture', scr, wave);
        
        
    end



    function tex = make_image(imagefile, colourmod, alpha)
        
        img = imread(imagefile);
        img = imresize(img, stimsize/size(img, 1));
        img = cat(3, rgb2gray(img), rgb2gray(img), rgb2gray(img));
        colourmod = colourmod/255;
        
        img(:,:,1) = img(:,:,1)*colourmod(1);
        img(:,:,2) = img(:,:,2)*colourmod(2);
        img(:,:,3) = img(:,:,3)*colourmod(3);
        
        
        if nargin < 3
            tex = Screen('MakeTexture', scr, img);
        elseif nargin == 3
            img(:,:,4) = alpha;
            tex = Screen('MakeTexture', scr, img);
        end
    
    end



    function take_break
        for k = 1:20
            Screen('DrawText', scr, 'Pause', xcen-offset-30, ycen-10, 255);
            Screen('DrawText', scr, 'Pause', xcen+offset-30, ycen-10, 255);
            Screen('DrawText', scr, num2str(k), xcen-5, 0, 255);
            Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
            Screen('Flip', scr);
            [~,~,pressList] = KbCheck;
            if pressList(ent_key)
                break
            elseif pressList(esc_key)
                error('You interrupted the script!');
            end
            WaitSecs(1);
        end
        
    end

    function take_long_break
            Screen('DrawText', scr, 'Long', xcen-offset-30, ycen-10, 255);
            Screen('DrawText', scr, 'Long', xcen+offset-30, ycen-10, 255);
            Screen('DrawText', scr, 'Pause', xcen-offset-30, ycen+10, 255);
            Screen('DrawText', scr, 'Pause', xcen+offset-30, ycen+10, 255);
            Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
            Screen('Flip', scr);

        
        for k = 1:60
            Screen('DrawText', scr, 'Long', xcen-offset-30, ycen-10, 255);
            Screen('DrawText', scr, 'Long', xcen+offset-30, ycen-10, 255);
            Screen('DrawText', scr, 'Pause', xcen-offset-30, ycen+10, 255);
            Screen('DrawText', scr, 'Pause', xcen+offset-30, ycen+10, 255);
            Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
            Screen('Flip', scr);
            [~,~,pressList] = KbCheck;
            if pressList(ent_key)
                break
            elseif pressList(esc_key)
                error('You interrupted the script!');
            end
            WaitSecs(1);
        end
        
    end


    function course = make_smooth_course(blockdur, type)
        events = load(['event_durations.mat']);
        course = zeros(1, ceil(blockdur*60));
        course(1) = randsample([1, 90], 1);
        
        if strcmp(type, 'gratings')
            events.adjDomA = events.adjDomA * 1.73/mean(events.adjDomA);
            events.adjDomC = events.adjDomC * 1.73/mean(events.adjDomC);
        end

        k_next = 1;
        for k = 1 : (ceil(blockdur*60))
            
        if k == k_next
            
            if rand > 0.2234 % proportion of reversions
                isSwitch = 1; % switch
            else
                isSwitch = 0; % reversion
            end
            
            
            while k_next == k && k < ceil(blockdur*60)-60
                mixDur = randsample([events.adjMixA; events.adjMixC], 1);
                domDur = randsample([events.adjDomA; events.adjDomC], 1);
                
                k_next = k + ceil(domDur*60) + ceil(mixDur*60) + mod(ceil(mixDur*60), 2);
            end
            
                        
            if k > ceil(blockdur*60)-60
                course(k_next : ceil(blockdur*60))=course(k_next);
            else
                course(k+1 : k+ceil(domDur*60)) = course(k)*ones(1, ceil(domDur*60));
                
                
                course(k+ceil(domDur*60)+1 : k+ceil(domDur*60)+ceil(mixDur*60)+mod(ceil(mixDur*60), 2)) = ...
                    round([linspace(course(k+ceil(domDur*60)), 45, ceil(mixDur*60/2)), linspace(45, abs(course(k)-isSwitch*90)+isSwitch, ceil(mixDur*60/2))]);
                
                
                
            end
            
        end
        end
        course(1 : find(course==45, 1)) = round(45 + 10*sin(2*pi* [1:find(course==45, 1)] / find(course==45, 1) ));
        
        course(ceil(blockdur*60) : size(course, 2)) = [];
    end


    function course = make_sudden_course(blockdur, type)
        events = load(['event_durations.mat']);
        course = zeros(1, ceil(blockdur*60));
        s = randsample([1, 90], 1);
        
        if strcmp(type, 'gratings')
            events.adjDomA = events.adjDomA * 1.73/mean(events.adjDomA);
            events.adjDomC = events.adjDomC * 1.73/mean(events.adjDomC);
        end
        
        k_next = 1;
        for k = 1 : (ceil(blockdur*60))
            
        if k == k_next
            
            if rand > 0.2234 % proportion of reversions
                isSwitch = 1; % switch
            else
                isSwitch = 0; % reversion
            end
            
            
            while k_next == k && k < ceil(blockdur*60)-60
                mixDur = randsample([events.adjMixA; events.adjMixC], 1);
                domDur = randsample([events.adjDomA; events.adjDomC], 1);
                
                k_next = k + ceil(domDur*60) + ceil(mixDur*60) + mod(ceil(mixDur*60), 2);
            end
            
                        
            if k > ceil(blockdur*60)-60
                course(k_next : ceil(blockdur*60))=course(k_next);
            else
                course(k+1 : k+ceil(domDur*60)) = s*ones(1, ceil(domDur*60));
                
                s = abs(s-isSwitch*90) + isSwitch*1;
                
                course(k+ceil(domDur*60)+1 : k+ceil(domDur*60)+ceil(mixDur*60)+mod(ceil(mixDur*60), 2)) = ...
                    45 * ones(1, ceil(mixDur*60)+mod(ceil(mixDur*60), 2));
                
                
                
            end
            
        end
        end
        course(1 : find(course==45, 1)) = 45;
        
        course(ceil(blockdur*60) : size(course, 2)) = [];
    end

%% Actual Experiment

    function [pressSecs, pressList] = rivalry_block(textures, blockdur, type)
        
        
        
        
        k = 0;
        pressSecs = zeros(blockdur/samplerate, 1);
        pressList = zeros(blockdur/samplerate, 3);
        
        if strcmp(type, 'images')
            
            Screen('FillRect', scr, [255 215 0]);
            Screen('FrameOval', scr, frameColour, frameRect, frameWidth);
            Screen('DrawText', scr, 'Red - Right', xcen - offset - stimsize/2, ycen - stimsize/4, 0);
            Screen('DrawText', scr, 'Green - Left', xcen - offset - stimsize/2, ycen + stimsize/5, 0);
            Screen('DrawText', scr, 'Red - Right', xcen + offset - stimsize/2, ycen - stimsize/4, 0);
            Screen('DrawText', scr, 'Green - Left', xcen + offset - stimsize/2, ycen + stimsize/5, 0);
            Screen('Flip', scr);
            WaitSecs(0.5);
            KbWait;
            Screen('FrameOval', scr, frameColour, frameRect, frameWidth);
            Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
            Screen('Flip', scr);
            trialOn = GetSecs + 2;
        
        
            PsychPortAudio('Start', pa, [], trialOn);
            Screen('DrawTexture', scr, textures(1), [], stimRect(1:4, 1));
            Screen('DrawTexture', scr, textures(2), [], stimRect(1:4, 2));
            Screen('FrameOval', scr, frameColour, frameRect, frameWidth);
            %Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
            Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
            t0 = Screen('Flip', scr, trialOn);
            last_wake = t0;

            while GetSecs < t0 + blockdur

                k = k+1;
                last_wake = WaitSecs('UntilTime', last_wake + samplerate);
                [~, pressSecs(k), firstPress] = KbCheck;

                pressList(k, 1:3) = firstPress(1, [l_key, u_key, r_key]);

                if firstPress(esc_key)
                    error('You interrupted the script!');
                end
            end

            pressSecs = pressSecs - t0;
            
        elseif strcmp(type, 'gratings')
            
            Screen('FillRect', scr, [127.5 127.5 127.5]);
            Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
            Screen('DrawText', scr, 'CW - Right', xcen - offset - stimsize/2, ycen - stimsize/4, 0);
            Screen('DrawText', scr, 'CCW - Left', xcen - offset - stimsize/2, ycen + stimsize/5, 0);
            Screen('DrawText', scr, 'CW - Right', xcen + offset - stimsize/2, ycen - stimsize/4, 0);
            Screen('DrawText', scr, 'CCW - Left', xcen + offset - stimsize/2, ycen + stimsize/5, 0);
            Screen('Flip', scr);
            WaitSecs(0.5);
            KbWait;
            Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
            Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
            Screen('Flip', scr);
            trialOn = GetSecs + 2;
        
        
            PsychPortAudio('Start', pa, [], trialOn);
            Screen('DrawTexture', scr, textures(1), [], gratRect(1:4, 1));
            Screen('DrawTexture', scr, textures(2), [], gratRect(1:4, 2));
            Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
            %Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
            Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
            t0 = Screen('Flip', scr, trialOn);
            last_wake = t0;

            while GetSecs < t0 + blockdur

                k = k+1;
                last_wake = WaitSecs('UntilTime', last_wake + samplerate);
                [~, pressSecs(k), firstPress] = KbCheck;

                pressList(k, 1:3) = firstPress(1, [l_key, u_key, r_key]);

                if firstPress(esc_key)
                    error('You interrupted the script!');
                end
            end

            pressSecs = pressSecs - t0;
            
        end
        
        
    end

    function [pressSecs, pressList, course] = rivalry_simulation(type, blockdur, LR, imgNrs, coursetype)
        
        load(['gaussians', num2str(randi(5)), '.mat']);
        
        if nargin > 4 && strcmp(coursetype, 'sudden')
            course = make_sudden_course(blockdur, type);
        else
            course = make_smooth_course(blockdur, type);
        end
        
        
        if strcmp(type, 'images');
            
            img1 = fullfile(pwd, 'Photos', [num2str(imgNrs(1)), '.jpg']);
            img2 = fullfile(pwd, 'Photos', [num2str(imgNrs(2)), '.jpg']);
            
            gaussians(1:2:180) = [];
            
            if LR == 1
                for i = 1:90
                    gaussians(i).red = imresize(gaussians(i).red, stimsize/size(gaussians(i).red, 1));
                    tex(i) = make_image(img1, [255 0 0], gaussians(i).red);
                end
                tex2 = make_image(img2, [0 255 0]);
            elseif LR == 2
                for i = 1:90
                    gaussians(i).red = imresize(gaussians(i).red, stimsize/size(gaussians(i).red, 1));
                    tex(i) = make_image(img1, [0 255 0], gaussians(i).red);
                end
                tex2 = make_image(img2, [255 0 0]);
            end
            
            
            Screen('FillRect', scr, [255 215 0]);
            Screen('FrameOval', scr, frameColour, frameRect, frameWidth);
            Screen('DrawText', scr, 'Red - Right', xcen - offset - stimsize/2, ycen - stimsize/4, 0);
            Screen('DrawText', scr, 'Green - Left', xcen - offset - stimsize/2, ycen + stimsize/5, 0);
            Screen('DrawText', scr, 'Red - Right', xcen + offset - stimsize/2, ycen - stimsize/4, 0);
            Screen('DrawText', scr, 'Green - Left', xcen + offset - stimsize/2, ycen + stimsize/5, 0);
            Screen('Flip', scr);
            WaitSecs(0.5);
            KbWait;
            
            Screen('FrameOval', scr, frameColour, frameRect, frameWidth);
            Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
            Screen('Flip', scr);
            trialOn = GetSecs + 2;
            PsychPortAudio('Start', pa, [], trialOn);
            Screen('DrawTexture', scr, tex(1), [], stimRect(1:4, 1));
            Screen('DrawTexture', scr, tex2, [], stimRect(1:4, 2));
            Screen('FrameOval', scr, frameColour, frameRect, frameWidth);
            Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
            t0 = Screen('Flip', scr, trialOn);
            last_flip = t0;
            
            
        elseif strcmp(type, 'gratings')
            
            gaussians(1:2:180) = [];
            
            if LR == 1
                for i = 1:90
                    gaussians(i).red = imresize(gaussians(i).red, stimsize/size(gaussians(i).red, 1));
                    tex(i) = make_grating(45, 0.8, gaussians(i).red);
                end
                tex2 = make_grating(-45, 0.8);
            elseif LR == 2
                for i = 1:90
                    gaussians(i).red = imresize(gaussians(i).red, stimsize/size(gaussians(i).red, 1));
                    tex(i) = make_grating(-45, 0.8, gaussians(i).red);
                end
                tex2 = make_grating(45, 0.8);
            end
            
            Screen('FillRect', scr, 127.5);
            Screen('FrameOval', scr, frameColour, frameRect, frameWidth);
            Screen('DrawText', scr, 'CW - Right', xcen - offset - stimsize/2, ycen - stimsize/4, 0);
            Screen('DrawText', scr, 'CCW - Left', xcen - offset - stimsize/2, ycen + stimsize/5, 0);
            Screen('DrawText', scr, 'CW - Right', xcen + offset - stimsize/2, ycen - stimsize/4, 0);
            Screen('DrawText', scr, 'CCW - Left', xcen + offset - stimsize/2, ycen + stimsize/5, 0);
            Screen('Flip', scr);
            WaitSecs(0.5);
            KbWait;
            
            Screen('FrameOval', scr, frameColour, frameRect, frameWidth);
            Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
            Screen('Flip', scr);
            trialOn = GetSecs + 2;
            PsychPortAudio('Start', pa, [], trialOn);
            Screen('DrawTexture', scr, tex(1), [], stimRect(1:4, 1));
            Screen('DrawTexture', scr, tex2, [], stimRect(1:4, 2));
            Screen('FrameOval', scr, frameColour, frameRect, frameWidth);
            Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
            t0 = Screen('Flip', scr, trialOn);
            last_flip = t0;
        end
            
            
            pressSecs = zeros(ceil(blockdur/samplerate), 1);
            pressList = zeros(ceil(blockdur/samplerate), 3);
            k = 0;
            l = 0;

            while GetSecs < t0 + blockdur
                
                l = l+1;
                Screen('DrawTextures', scr, [tex2, tex2, tex(course(l)), tex(course(l))], [], [stimRect(1:4, 1:2), stimRect(1:4, 1:2)]);
                if strcmp(type, 'images')
                    Screen('FrameOval', scr, frameColour, frameRect, frameWidth);
                elseif strcmp(type, 'gratings')
                    Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
                end
                Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
                last_flip = Screen('Flip', scr, last_flip + 0.015);
                if l == 1;
                imArray = Screen('GetImage', scr, stimRect(1:4, 1));
                imwrite(imArray, [ type, '.png' ]);
                end
                
                for i = 0:3
                    k = k+1;
                    WaitSecs('UntilTime', last_flip + i*samplerate);
                    [~, pressSecs(k), firstPress] = KbCheck;
                    pressList(k, 1:3) = firstPress(1, [l_key, u_key, r_key]);
                    if firstPress(esc_key)
                        error('You interrupted the script!');
                    end
                end
                
            end
            
                if strcmp(type, 'images')
                    Screen('FrameOval', scr, frameColour, frameRect, frameWidth);
                elseif strcmp(type, 'gratings')
                    Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
                end
                Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
                last_flip = Screen('Flip', scr);
            
            pressSecs = pressSecs - t0;
            Screen('Close', tex);
            Screen('Close', tex2);
    end

    function demonstrate_stimuli(textures, type)
        
        if strcmp(type, 'images');
            
        Screen('FillRect', scr, [255 215 0]);
        Screen('Flip', scr);
        
        Screen('DrawTexture', scr, textures(1), [], stimRect(1:4, 1));
        Screen('DrawTexture', scr, textures(1), [], stimRect(1:4, 2));
        Screen('FrameOval', scr, frameColour, frameRect, frameWidth);
        Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
        Screen('Flip', scr);
        
        WaitSecs(1); KbWait;
        
        Screen('DrawTexture', scr, textures(2), [], stimRect(1:4, 1));
        Screen('DrawTexture', scr, textures(2), [], stimRect(1:4, 2));
        Screen('FrameOval', scr, frameColour, frameRect, frameWidth);
        Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
        Screen('Flip', scr);
        WaitSecs(1); KbWait;
        
        Screen('DrawTexture', scr, textures(1), [], stimRect(1:4, 1));
        Screen('DrawTexture', scr, textures(2), [], stimRect(1:4, 2));
        Screen('FrameOval', scr, frameColour, frameRect, frameWidth);
        Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
        Screen('Flip', scr);
        imArray = Screen('GetImage', scr);
        imwrite(imArray, 'object_stimuli.png');
        
        WaitSecs(1); KbWait;
        
        Screen('FrameOval', scr, frameColour, frameRect, frameWidth);
        Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
        Screen('Flip', scr);
        
        WaitSecs(1); KbWait;
        
        elseif strcmp(type, 'gratings');
        
        Screen('FillRect', scr, 127.5);
        Screen('Flip', scr);
        
        Screen('DrawTexture', scr, textures(1), [], gratRect(1:4, 1));
        Screen('DrawTexture', scr, textures(1), [], gratRect(1:4, 2));
        Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
        Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
        Screen('Flip', scr);
        
        WaitSecs(1); KbWait;
        
        Screen('DrawTexture', scr, textures(2), [], gratRect(1:4, 1));
        Screen('DrawTexture', scr, textures(2), [], gratRect(1:4, 2));
        Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
        Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
        Screen('Flip', scr);
        WaitSecs(1); KbWait;
        
        Screen('DrawTexture', scr, textures(1), [], gratRect(1:4, 1));
        Screen('DrawTexture', scr, textures(2), [], gratRect(1:4, 2));
        Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
        Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
        Screen('Flip', scr);
        imArray = Screen('GetImage', scr);
        imwrite(imArray, 'grating_stimuli.png');
        
        WaitSecs(1); KbWait;
        
        Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
        Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
        Screen('Flip', scr);
        WaitSecs(1); KbWait;
        end
        
    end


end
