%% Manual line bisection task - AG.Mitchell 18.08.18

clc
clear all
datenow(1:6) = fix(clock);
dummymode = 0;
practice = 0 ;
useTouch = true;

%% Variables
matfilename = sprintf('MLB_%d%d%d_%d%d%d', datenow);
% Struct
stim = struct;
stim.line = [];
stim.time = [];
stim.bisect = [];
response = struct;
response.touch = zeros(180,2); %predefining response click
time = struct;
data = struct;

% Dimensions in cm
sd = 57; %distance from laptop
sx = 47.5; %size of screen
sy = 26.5;

% Size of stimuli
stim.line.width = 0.9; %line width in degrees
stimlinewdth_cm = tan(pi*stim.line.width/360)*2*sd; %in cm
stim.line.offset = 2; %offset left or right by 1.5 degrees
stimlineoff_cm = tan(pi*stim.line.offset/360)*2*sd;
% Different linelengths
stim.line.lengths = [10, 20, 30];
stimlinelength_cm = tan(pi*stim.line.lengths/360)*2*sd;  

% Import data matrix
[Data,Text] = xlsread('TrialMatrix_MLB'); %importing trial information for landmark tasks
stim.side = Data(:, ismember(Text, 'side')); %5deg left, 0, 5deg right
stim.size = Data(:, ismember(Text, 'size'));
stim.mouse = Data(:, ismember(Text, 'mouse')); %starting location of mouse, right or left of line

% Deciding which side of line is white/black to place  shift
upper = zeros(90,1); upper = upper+1; lower = upper+1;
position = [upper; lower];
stim.pos = position(randperm(length(position))); %here this is used to randomise contrast order

if practice == 0
    nrtrials = length(stim.side); %getting number of trials
else
    nrtrials = 10;
end

% Timing
stim.time.response = 5; %length of response in s
stim.time.fix = 1; %fixation 500ms
stim.time.line = 0.2; %length of line onset

%% Output
response = struct;
response.actual = [];
response.expect = [];
        
%% Starting 
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
inc = white-grey;

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

% Making stimuli
stim.line.lengthspix = stimlinelength_cm*vadxcm; %pixels, 
stim.line.widthpix = stimlinewdth_cm*vadycm; 
%stim.dim.bisectlength.pix = stim.dim.bisectlength.cm*DisplayYSize/sy; %verticle so along y axis
stim.line.offsetpix = stimlineoff_cm*vadxcm; 
%stim.dim.dist.pix = stim.dim.dist.cm*DisplayXSize/sx;

%% Get touchscreen info
if useTouch
    touch.device = GetTouchDeviceIndices(); %getting device for touchscreen
    slots = 1; %number of 'touches' at one point
    touch.waitsecs = 5 ; %maximum response time before moving on
    % Create and start touch queue for window and touchscreen:
    TouchQueueCreate(window, touch.device, [], [], [], 8 );
    TouchQueueStart(touch.device); 
end

%% Starting experiment
% Welcome screen
Screen('FillRect', window, grey);
Screen('TextStyle', window, 0); %normal
Screen('TextFont', window, 'Lucida Console');
text = 'Press a button to start bisection task';
width = RectWidth(Screen ('TextBounds', window, text));
Screen('DrawText', window, text, DisplayXSize/2 - width/2, DisplayYSize/2, black); 
Screen('Flip', window); 

while KbCheck
end %wait for all keys are released, then wait for the first key press
keyisdown = 0;
while ~keyisdown
    [keyisdown] = KbCheck;
    WaitSecs(0.001); % delay to prevent CPU hogging
end

text = 'Touch the middle of the line';
width = RectWidth(Screen ('TextBounds', window, text));
Screen('DrawText', window, text, DisplayXSize/2 - width/2, DisplayYSize/2, black); 
Screen('Flip', window); 
WaitSecs(5)
Screen('Flip', window);

% Loop starts
time.experiment.start = GetSecs;
for i = 1:nrtrials  
    % Stimulus line lengths
    switch stim.size(i)
        case 1
            linelength = stim.line.lengthspix(1);
        case 2
            linelength = stim.line.lengthspix(2);
        case 3
            linelength = stim.line.lengthspix(3);
    end
    
    % switching between sides of shift
    switch stim.side(i)
        case 0
            offset = 0;
        case 1
            offset = -stim.line.offsetpix;
        case 2
            offset = stim.line.offsetpix;
    end
    
    % Coordinates for horizontal line
    liney = midY; %position along y axis
    linex1 = midX-linelength/2;
    linex2 = midX+linelength/2;
    
    % deciding where the mouse will start before each trial    
    switch stim.mouse(i)
        case 1
            mstart = round(linex1);
        case 2
            mstart = round(linex2);
    end
    
    %% Present stim
    text = '+';
    width = RectWidth(Screen ('TextBounds', window, text));
    height = RectHeight(Screen ('TextBounds', window, text));
    Screen('Flip', window);
    Screen('DrawText', window, text, midX, midY-height/2, black); %fixation cross at center
    Screen('Flip', window);
    WaitSecs(stim.time.fix)
    

   % if stim.side(i) == 1 %if side longer on left
        % Draw one line with different offset
        if stim.pos(i) == 1 %top: black, white; bottom: white, black
            Screen('FillRect', window, white, [linex1+offset, liney-(stim.line.widthpix/2), linex2+offset, liney]); %line 1 - top left
            Screen('FillRect', window, black, [linex1+offset, liney, linex2+offset, liney+(stim.line.widthpix/2)]); % line 2 - bottom right
        elseif stim.pos(i) == 2 %top: white, black; bottom: black, white
            Screen('FillRect', window, black, [linex1+offset, liney-(stim.line.widthpix/2), linex2+offset, liney]); %line 1 - top left
            Screen('FillRect', window, white, [linex1+offset, liney, linex2+offset, liney+(stim.line.widthpix/2)]); % line 2 - bottom right
        end
   % end
   
   % Create and start touch queue for window and touchscreen:


    Screen('Flip', window);
    time.target.start(i) = GetSecs;
    %WaitSecs(stim.time.line)
    if useTouch
        TouchEventFlush(touch.device)
        event = TouchEventGet(touch.device, window, touch.waitsecs);
        navail = TouchEventAvail(touch.device);
    else
        WaitSecs(2)
    end
    %end
    
    %KbStrokeWait; %for testing
    Screen('Flip',window);
    time.target.end(i) = GetSecs;
    time0=GetSecs; %for response
    %KbStrokeWait; 
    
    % Getting expected response
    if isstruct(event) == 1 %if receive an event1
        response.touch(i,1) = event.MappedX;
        response.touch(i,2) = event.MappedY;       
    else
        response.touch(i,1) = NaN;
        response.touch(i,2) = NaN; 
    end
    % Getting actual line middle
    response.middle(i,:) = midX+offset;
    % Releasing queue
    
    if i == (length(nrtrials)/2);
        text = 'Switch hands...';
        width = RectWidth(Screen ('TextBounds', window, text));
        Screen('DrawText', window, text, DisplayXSize/2 - width/2, DisplayYSize/2, black); 
        Screen('Flip', window); 
        KbStrokeWait;
        Screen('Flip', window);
    end
end
time.experiment.end = GetSecs;
TouchQueueRelease(touch.device)

% Thank you screen
    thankYouText = 'End of task, thank you!';
    height=RectHeight(Screen ('TextBounds', window, thankYouText));
    width=RectWidth(Screen ('TextBounds', window, thankYouText));
    Screen('DrawText', window, thankYouText, DisplayXSize/2 - width/2, DisplayYSize/2, black);
    Screen('Flip', window);
    WaitSecs (2);
    Screen('Flip', window);
    Screen('CloseAll');

%% Analysing
% Results
 %data.RT = response.RT*1000; %transforming to ms
load('ts_calib.mat');
% correcting response with calibration results
response.touchcorr(:,1) = response.touch(:,1) - calResponse.avError(1);
response.touchcorr(:,2) = response.touch(:,2) - calResponse.avError(2);
% 
time.target.dur = time.target.end - time.target.start;
time.experiment.dur = time.experiment.end - time.experiment.start;
% 
% %% Combining results for analysis
response.middle = response.middle/vadx; %converting results from pix to degrees
response.touchcorr = response.touchcorr/vadx;
response.error = response.touchcorr-response.middle; %everything -ve left-hand error, everything +ve right hand error

% removing outlier trials from the experiment
averror = nanmean(response.error);
sderror = nanstd(response.error);
response.error(find(response.error > (averror+2.5*sderror))) = NaN;
response.error(find(response.error < (averror-2.5*sderror))) = NaN;
data.matrix = [stim.size, stim.side, response.middle, response.touch, response.error]; %matrix of results 
% removing any NaN trials (2.5 sds away from mean and no response trials)
data.matrix = data.matrix(find(data.matrix(:,5) >= 0.0001|data.matrix(:,5) <= 0.0001),:); %removing NaN trials

% remaining mean and std
data.error.mean = nanmean(response.error);
data.error.sd = nanstd(response.error);

%% Don't forget to get read of NaN values within matrix
% 
% data.left = data.matrix(find(data.matrix(:,7) <0),:);
% data.right = data.matrix(find(data.matrix(:,7) >0),:);

%% Saving data
catch
    rethrow(lasterror) 
    sca
    save(matfilename, 'stim', 'response', 'time', 'data')
    TouchQueueRelease(touch.device)
end
sca
save(matfilename, 'stim', 'response', 'time', 'data')