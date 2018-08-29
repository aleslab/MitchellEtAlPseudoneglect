%% Touch screen calibration
%% AG. Mitchell 27.08.18
% This code runs a calibration of the touchscreen manually
% Presents a cross at centre folowed by 12 points on the screen, waits for a 'touch' input.
% Cross changes colour once accurately touched, then disappears and the
% next cross is presented. Error is measured (up, down, left or right) for
% each calibration point, if it's over 0.5deg then calibration is redone. 
% Centre cross presented at beginning and end of calibration

clear all
clc
datenow(1:6) = fix(clock);
dummymode = 1;
useTouch = true;

%% Variables

% Filename
matfilename = sprintf('touchscreen_cal_%d%d%d%d%d%d.mat', datenow);

% save data
calStim = struct;
calResponse = struct;
calStim.text = '+';

% Dimensions in cm
sd = 57; %distance from laptop
sx = 47.5; %size of screen
sy = 26.5;% Dimensions in cm

%% Setting up ptb
try
    % Variables
    screens = Screen('Screens');
    screenNumber = max(screens);
    %Basically gets rid of all sync warnings if not needed
    Screen('Preference', 'SkipSyncTests', 1);
    Screen('Preference','SuppressAllWarnings', 1);
    Screen('Preference','VisualDebugLevel', 0);

    %define back and white (white = 1, black = 0)
    white = WhiteIndex(screenNumber);
    black = BlackIndex(screenNumber);
    grey=white/2; %white is too strong on CRT
    red = [1 0 0];

    if dummymode==1%For debugging use a smaller window
        [window,rect] = Screen('OpenWindow', screenNumber, grey,[RectLeft, RectTop RectLeft+850 RectTop+250]); %for testing    
    else
        [window,rect]= Screen('OpenWindow', screenNumber, grey); %full screen window
    end
    [winWidth, winHeight] = WindowSize(window); 
    
    % Text size
    textpt = 0.0352; %cm for 1 pt
    textsizedeg = 1; %degrees, visual crowding in fovea under 0.2 deg; 0.35 deg visual angle, size of one letter according to Spinelli 2002;
    textsizecm = tan(pi*textsizedeg/360)*2*sd; %in cm
    textsizept = (round(textsizecm/textpt)); %in pt, for some reason this number needs to be doubled
    Screen('TextSize', window, textsizept);

    % Display dimensions
    DisplayXSize = rect(3);
    DisplayYSize = rect(4);
    midX = DisplayXSize/2;
    midY = DisplayYSize/2; 
    %find pixels/cm
    vadxcm=DisplayXSize/sx; %pix per 1 cm
    vadycm=DisplayYSize/sy;
    %find pixels/degree
    vadx=vadxcm/(atan(1/sd)*180/pi); % pixels for 1 degree visual angle
    vady=vadycm/(atan(1/sd)*180/pi); % pixels for 1 degree visual angle
    
    %% Get touchscreen info
    if useTouch
        touch.device = GetTouchDeviceIndices(); %getting device for touchscreen
        slots = 1; %number of 'touches' at one point
        touch.waitsecs = 10; %maximum response time before moving on
    end

    %% Setting calibration points
    % Making matrix of calibration coordinates to run through (relative to
    % centre)
    width = RectWidth(Screen ('TextBounds', window, calStim.text));
    height = RectHeight(Screen ('TextBounds', window, calStim.text));
    calStim.mat(1,:) = [midX, midY]; %centre
    calStim.mat(2,:) = [0+width, 0+height]; %top left coords
    calStim.mat(3,:) = [midX, 0+height]; %top middle
    calStim.mat(4,:) = [DisplayXSize-width, 0+height]; %top right
    calStim.mat(5,:) = [DisplayXSize-width, midY]; %cent right
    calStim.mat(6,:) = [DisplayXSize-width, DisplayYSize-height]; %bottom right
    calStim.mat(7,:) = [midX, DisplayYSize-height]; %bottom middle
    calStim.mat(8,:) = [0+width, DisplayYSize-height]; %bottom left
    calStim.mat(9,:) = [0+width, midY]; %cent left
    calStim.mat(10,:) = [midX/2, midY/2]; %small screen top left
    calStim.mat(11,:) = [midX + midX/2, midY/2]; %small screen top right
    calStim.mat(12,:) = [midX + midX/2, midY + midY/2]; %small screen bottom right
    calStim.mat(13,:) = [midX/2, midY + midY/2]; %small screen bottom left
    calStim.mat(14,:) = [midX, midY]; %centre
    
    %% Starting calibration
    line = 'Touchscreen calibration, press any key to start';
    widthl = RectWidth(Screen ('TextBounds', window, line));
    Screen('DrawText', window, line, DisplayXSize/2 - widthl/2, DisplayYSize/2, black); 
    Screen('Flip', window); 
    KbStrokeWait;
    Screen('Flip', window);
    
    %% Running calibration
    for i = 1:length(calStim.mat)
        Screen('DrawText', window, calStim.text, calStim.mat(i,1), calStim.mat(i,2), black);
        Screen('Flip', window)

        if useTouch
            % Create and start touch queue for window and touchscreen:
            TouchQueueCreate(window, touch.device);
            TouchQueueStart(touch.device);
            
            event = TouchEventGet(touch.device, window, touch.waitsecs);
        else
            WaitSecs(2)
        end
        
                %% Getting response
        if useTouch
            if isstruct(event) %incase there is no response
                calResponse.actual(i,1) = event.MappedX; 
                calResponse.actual(i,2) = event.MappedY;
            else
                calResponse.actual(i,1) = NaN; 
                calResponse.actual(i,2) = NaN;
            end
        else           
            calResponse.actual(i,1) = NaN; 
            calResponse.actual(i,2) = NaN;
        end
        
        errorX(i) = calResponse.actual(i,1) - calStim.mat(i,1); %square rooted so that the 
        errorY(i) = calResponse.actual(i,2) - calStim.mat(i,2);
     
        Screen('DrawText', window, calStim.text, calStim.mat(i,1), calStim.mat(i,2), red);
        Screen('Flip', window);
        WaitSecs(0.6)
                 
    end
    sca
catch
    rethrow(lasterror)
    sca
    save(matfilename)
end
save(matfilename)
