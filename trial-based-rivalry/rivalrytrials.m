
function rivalrytrials
try
clear all
%% Session Variables
commandwindow;
ID = input('Participant ID? ', 's');
scr_diagonal = 24; % in inches
scr_distance = 60; % in cm
% diagnosis = input('Diagnosis? ');

tstamp = clock;
if ~isdir( fullfile(pwd, 'Results', 'rivalry trials repetition') )
    mkdir( fullfile(pwd, 'Results', 'rivalry trials repetition') );
end
savefile = fullfile(pwd, 'Results', 'rivalry trials repetition', [sprintf('%02d-%02d-%02d-%02d%02d-', tstamp(1), tstamp(2), tstamp(3), tstamp(4), tstamp(5)), ID, '.mat']);


%% Experiment Variables

scr_background = 127.5;
scr_no = max(Screen('Screens'));
scr_dimensions = Screen('Rect', scr_no);
xcen = scr_dimensions(3)/2;
ycen = scr_dimensions(4)/2;

stimsizedeg = 3.5;
stimsize = visual_angle2pixel(stimsizedeg, scr_diagonal, scr_distance, scr_no);

trialdur = 6;
samplerate = 0.005; % time in between sample collection, in seconds

iti = 1; % MINIMUM time in between trials

reps = 9; % repitions per condition - this means per counterbalanced rep of Left/Right swap & Green/Red Exposure



%% Set-Up

% Keyboard
KbName('UnifyKeyNames');
l_key = KbName('LeftArrow'); r_key = KbName('RightArrow');
u_key = KbName('UpArrow'); d_key = KbName('DownArrow');
esc_key = KbName('Escape');
ent_key = KbName('Return'); ent_key = ent_key(1);

% Sound
InitializePsychSound;
pa = PsychPortAudio('Open', [], [], [], [], [], 256);
bp400 = PsychPortAudio('CreateBuffer', pa, [MakeBeep(400, 0.2); MakeBeep(400, 0.2)]);
PsychPortAudio('FillBuffer', pa, bp400);

% Screen
scr = Screen('OpenWindow', scr_no, scr_background);
HideCursor;
if IsLinux % Correct font for Linux for some reason...
    Screen('TextFont', scr, '-schumacher-clean-medium-r-normal--0-0-75-75-c-0-koi8-r');
end
Screen('BlendFunction', scr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Save script with data for posterity
scripts = savescripts;

% Priority - ensure this runs with best timing
Priority(1);



%% Stimuli

% Frame
frameWidth = 5;
frameColour = 0;

% Fixation
fixWidth = stimsize/30 + mod(stimsize/30, 2)-1;
fixLength = stimsize/15 + mod(stimsize/15, 2);
fixColour = 255;



%% Find merge spot
% In this routine, increase proximity (by pressing UP) of two boxes until the participant
% reports that the "edges meet", i.e. the boxes are immediately adjacent.
% Then press ENTER to make the boxes jump, match color - this will make
% them merge for participant. Colored letters are inside the merged box for
% confirmation (ask if they can see letters fine). Then press ENTER to
% accept the calibration and move on.
[offset, frameRect, stimRect, fixLines] = find_offset;
Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
Screen('Flip', scr);
WaitSecs(0.5);
KbStrokeWait;



%% Demonstrate Stimuli
textures(1) = make_image(fullfile(pwd, 'Photos', '1.jpg'), [255 0 0]);
textures(2) = make_image(fullfile(pwd, 'Photos', '2.jpg'), [0 255 0]);

% This subfunction presents: Image 1 to both eyes, then Image 2 to both
% eyes, then Image 1&2 to both eyes. Talk Participant through this.
demonstrate_stimuli(textures);

Screen('Close', textures);


%% A block practice
[prac.LR, prac.exposeLR, prac.exposeDur, prac.exposeIsi, prac.catchTrial, prac.exposeR] = make_trial_schedule(2.0, 0.2, 1);
prac.n = size(prac.LR, 2);

for j = 1:size(prac.LR, 2)
    prac.imgNrs(j, 1:2) = randsample(15, 2);
    
    image1 = fullfile(pwd, 'Photos', [num2str(prac.imgNrs(j, 1)), '.jpg']);
    image2 = fullfile(pwd, 'Photos', [num2str(prac.imgNrs(j, 2)), '.jpg']);
    
    textures(prac.LR(j)) = make_image(image1, [255, 0, 0]);
    textures(mod(prac.LR(j), 2)+1) = make_image(image2, [0, 255, 0]);
    
    WaitSecs(iti);
    
        [prac.pressSecs(j).list(:,:), prac.pressList(j).list(:,:)] = trial_rivalry(textures, trialdur, prac.exposeLR(j), prac.exposeDur(j), prac.exposeIsi(j), prac.catchTrial(j));
        
        if prac.pressSecs(j).list(1,1) == -1
            break
        end
    
    Screen('Close', textures);
    
end

take_break;



%% Trial Block Pictures
    [data.LR, data.exposeLR, data.exposeDur, data.exposeIsi, data.catchTrial, data.exposeIm] = make_trial_schedule(2.0, 0.2, reps);
    data.n = size(data.LR, 2);
    
for j = 1:size(data.LR, 2)
    
    data.imgNrs(j, 1:2) = randsample(15, 2)';
    
    image1 = fullfile(pwd, 'Photos', [num2str(data.imgNrs(j, 1)), '.jpg']);
    image2 = fullfile(pwd, 'Photos', [num2str(data.imgNrs(j, 2)), '.jpg']);
    
    textures(data.LR(j)) = make_image(image1, [255, 0, 0]);
    textures(mod(data.LR(j), 2)+1) = make_image(image2, [0, 255, 0]);
            
    WaitSecs(iti);
    
    [data.pressSecs(j).list(:,:), data.pressList(j).list(:,:)] = trial_rivalry(textures, trialdur, data.exposeLR(j), data.exposeDur(j), data.exposeIsi(j), data.catchTrial(j));
    
    if data.pressSecs(j).list(1,1) == -1
        break
    end
    Screen('Close', textures);
    
    if j < size(data.LR, 2) && mod(j, 16) == 0
        take_break;
    end
end




%% Shutdown
Screen('DrawText', scr, 'Thank you!', xcen-offset-30, ycen-10, 255);
Screen('DrawText', scr, 'Thank you!', xcen+offset-30, ycen-10, 255);
Screen('Flip', scr);
WaitSecs(2);


% some quick post-block analysis of the trials.
for j = 1:data.n
    data.pressAv(:,:,j) = mean(data.pressList(j).list, 1);
    h = 0;
    data.switches(j) = 0;
    data.reversions(j) = 0;
    for jj = 2:size(data.pressSecs(j).list, 1)
        if   ~(isequal(data.pressList(j).list(jj, 1:3), [0 0 0]) || sum(data.pressList(j).list(jj, 1:3))>1) && ~isequal(data.pressList(j).list(jj, 1:3), data.pressList(j).list(jj-1, 1:3))
            h = h+1;
            data.transitionPress(h, 1:3, j) = data.pressList(j).list(jj, 1:3);
            data.transitionSecs(h, 1, j) = data.pressSecs(j).list(jj, 1);
        end
    end
    
    try
        start_ind = min([find(data.transitionPress(:, 1, j), 1), find(data.transitionPress(:, 3, j), 1)]);
    catch
        break
    end
    
    for jj = 3:max(find(data.transitionSecs(:, :, j), 1, 'last'))
        
        if isequal(data.transitionPress(jj-2:jj, 1:3, j), [1 0 0; 0 1 0; 0 0 1]) || isequal(data.transitionPress(jj-1:jj, 1:3, j), [1 0 0; 0 0 1])
            data.switches(j) = data.switches(j) + 1;
            data.domDur1(ceil((data.switches(j))/2), j) = data.transitionSecs(jj, 1, j) - data.transitionSecs(start_ind, 1, j);
            start_ind = jj;
        end
        
        if isequal(data.transitionPress(jj-2:jj, 1:3, j), [0 0 1; 0 1 0; 1 0 0]) || isequal(data.transitionPress(jj-1:jj, 1:3, j), [0 0 1; 1 0 0])
            data.switches(j) = data.switches(j) + 1;
            data.domDur3(ceil((data.switches(j)/2)), j) = data.transitionSecs(jj, 1, j) - data.transitionSecs(start_ind, 1, j);
            start_ind = jj;
        end
        
        if isequal(data.transitionPress(jj-2:jj, 1:3, j), [1 0 0; 0 1 0; 1 0 0]) || isequal(data.transitionPress(jj-2:jj, 1:3, j), [0 0 1; 0 1 0; 0 0 1])
            data.reversions(j) = data.reversions(j) + 1;
        end
    end
    
    if find(data.transitionPress(:, 1, j), 1, 'last') > find(data.transitionPress(:, 3, j), 1, 'last')
        data.domDur1(ceil((data.switches(j)+1)/2), j) = data.transitionSecs(find(data.transitionPress(:, 1, j), 1, 'last')) - data.transitionSecs(start_ind);
    elseif find(data.transitionPress(:, 3, j), 1, 'last') > find(data.transitionPress(:, 1, j), 1, 'last')
        data.domDur3(ceil((data.switches(j)+1)/2), j) = data.transitionSecs(find(data.transitionPress(:, 3, j), 1, 'last')) - data.transitionSecs(start_ind);
    end
    
end
data.pressAv = permute(data.pressAv, [2, 3, 1]);



% Close all windows, save data, return to normal matlab
sca;
save(savefile);
ListenChar(0);
Priority(0);



catch err
%% Catch
% Close everything, save data with Error marker, return to normal matlab
sca;
savefile = [savefile(1:(size(savefile, 2)-4)), '-ERROR.mat'];
save(savefile);
ListenChar(0);
Priority(0);
rethrow(err);



end
%% Subfunctions

%% Calibration, Break and Stimulus Creation

    function [offset, frameRect, stimRect, fixLines] = find_offset
        offset = stimsize;
        % Firstly, present white & black circles to avoid merging and find
        % the spot where they meet!
        while 1
            frameRect= [-stimsize/2-frameWidth;
                       -stimsize/2-frameWidth;
                        stimsize/2+frameWidth;
                        stimsize/2+frameWidth];
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
            [~, keyCode] = KbStrokeWait;
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
                error('You interrupted the script');
            end
            if offset < stimsize
                offset = stimsize;
            end
        end
        % Once the point where they meet has been established, make a jump
        % and then present the circles the same colour. Also present text
        % to verify
        while 1
            frameRect= [-stimsize/2-frameWidth;
                       -stimsize/2-frameWidth;
                        stimsize/2+frameWidth;
                        stimsize/2+frameWidth];
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
            [~, keyCode] = KbStrokeWait;
            
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
        
    end




    function tex = make_image(imagefile, colourmod)
        
        img = imread(imagefile);
        img = cat(3, rgb2gray(img), rgb2gray(img), rgb2gray(img));
        colourmod = colourmod/255;
        
        img(:,:,1) = img(:,:,1)*colourmod(1);
        img(:,:,2) = img(:,:,2)*colourmod(2);
        img(:,:,3) = img(:,:,3)*colourmod(3);
        
        tex = Screen('MakeTexture', scr, img);
    end

    

    function take_break
        for k = 1:30
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
        Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
        Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
        Screen('Flip', scr);
        KbStrokeWait;
    end


%% Actual Experiment


    function [pressSecs, pressList] = trial_rivalry(textures, duration, preExpose, preDuration, isi, catchTrial)
        
        k = 0;
        pressSecs = zeros(ceil(duration/samplerate), 1);
        pressList = zeros(ceil(duration/samplerate), 3);
        trialOn = GetSecs + 1.5;
        started = 0;
        
        % Exposure
        if nargin > 2 && preExpose
            
            if ~catchTrial
                % Texture for adaptation
                Screen('DrawTexture', scr, textures(preExpose), [], stimRect(1:4, preExpose));
                Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
                Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
                preOn = Screen('Flip', scr);
            elseif catchTrial
                % Present the tex for adaptation to the OPPOSITE eye
                Screen('DrawTexture', scr, textures(preExpose), [], stimRect(1:4, mod(preExpose, 2) + 1));
                Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
                Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
                preOn = Screen('Flip', scr);
            end
            
            % empty frame (ISI)
            Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
            Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
            preOff = Screen('Flip', scr, preOn+preDuration-1/120);
            trialOn = preOff + isi;
            
        end
                
        % both objects, with a beep
        PsychPortAudio('Start', pa, [], trialOn);
        Screen('DrawTexture', scr, textures(1), [], stimRect(1:4, 1));
        Screen('DrawTexture', scr, textures(2), [], stimRect(1:4, 2));
        Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
        Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
        t0 = Screen('Flip', scr, trialOn-1/120);
        last_wake = t0;
        tEnd = t0 + duration;
        
        while GetSecs < tEnd
            k = k + 1;
            
            WaitSecs('UntilTime', last_wake + samplerate);
            [~, pressSecs(k), singlePress] = KbCheck;
            last_wake = pressSecs(k);
            
            pressList(k, 1:3) = singlePress(1, [l_key, u_key, r_key]);
            
            if singlePress(esc_key)
                % Break out of script
                error('You interrupted the script!');
            elseif singlePress(ent_key)
                % Break out of loop
                pressList = zeros(duration/samplerate, 1)-1;
                pressSecs = zeros(duration/samplerate, 1)-1;
                return
            end
                
                if ~started && (pressList(k, 1) || pressList(k, 3))
                    started = 1;
                    tEnd = GetSecs + duration;
                end
                
            
        end
        
        % Return to just frame only
        Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
        Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
        Screen('Flip', scr);
        
        pressSecs = pressSecs - t0;
        % if nargin > 2 && preExpose && abs(t0-(trialOn+1/120)) > 1/120 || abs(preOff-(preOn+preDuration+1/120)) > 1/120
            % warning('Demanded and reported timing were off: t0-trialOn=%d, preOff-PreOn-preDuration=%d', t0-(trialOn+1/120), preOff-(preOn+preDuration+1/120));
        % end
    end


    function [LR, exposeLR, exposeDur, exposeIsi, catchTrial, exposeR] = make_trial_schedule(exposedurList, isiList, trialN)
        
        k = 0;
        % real trials
        for ilr = 1:2
            for iexposelr = 1:2
                for iexposedur = exposedurList
                    for iisi = isiList
                        for itrial = 1:trialN
                            k = k+1;
                            LR(k) = ilr;
                            exposeLR(k) = iexposelr;
                            exposeDur(k) = iexposedur;
                            exposeIsi(k) = iisi;
                        end
                    end
                end
            end
        end
        % catch trials
        for ilr = 1:2
            for iexposelr = 1:2
                for iexposedur = max(exposedurList)
                    for iisi = min(isiList)
                        for itrial = 1:trialN
                            k = k+1;
                            LR(k) = ilr;
                            exposeLR(k) = iexposelr;
                            exposeDur(k) = iexposedur;
                            exposeIsi(k) = iisi;
                            catchTrial(k) = 1;
                        end
                    end
                end
            end
        end
        % zero exposure trials
        for ilr = 1:2
            for iexposelr = 1:2
                for iexposedur = max(exposedurList)
                    for iisi = min(isiList)
                        for itrial = 1:ceil(trialN/2)
                            k = k+1;
                            LR(k) = ilr;
                            exposeLR(k) = 0;
                            exposeDur(k) = iexposedur;
                            exposeIsi(k) = iisi;
                            catchTrial(k) = 2;
                        end
                    end
                end
            end
        end
        
        sched = randperm(k);
        LR = LR(sched);
        exposeLR = exposeLR(sched);
        exposeDur = exposeDur(sched);
        exposeIsi = exposeIsi(sched);
        catchTrial = catchTrial(sched);
        exposeR = (LR == exposeLR);
        
    end

    
    function demonstrate_stimuli(textures)
        Screen('DrawTexture', scr, textures(1), [], stimRect(1:4, 1));
        Screen('DrawTexture', scr, textures(1), [], stimRect(1:4, 2));
        Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
        Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
        Screen('Flip', scr);
        KbStrokeWait;
        
        Screen('DrawTexture', scr, textures(2), [], stimRect(1:4, 1));
        Screen('DrawTexture', scr, textures(2), [], stimRect(1:4, 2));
        Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
        Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
        Screen('Flip', scr);
        KbStrokeWait;
        
        Screen('DrawTexture', scr, textures(1), [], stimRect(1:4, 1));
        Screen('DrawTexture', scr, textures(2), [], stimRect(1:4, 2));
        Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
        Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
        Screen('Flip', scr);
        KbStrokeWait;
        
        Screen('FrameRect', scr, frameColour, frameRect, frameWidth);
        Screen('DrawLines', scr, fixLines, fixWidth, fixColour);
        Screen('Flip', scr);
        KbStrokeWait;
        
    end


end
