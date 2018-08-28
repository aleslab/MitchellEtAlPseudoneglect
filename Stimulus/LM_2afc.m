%% Landmarks task - AG.Mitchell 09.05.18

clc
clear all
datenow(1:6) = fix(clock);
dummymode = 0;
practice = 0;

%% Variables
matfilename = sprintf('LM_%d%d%d_%d%d%d', datenow);
% Struct
stim = struct;
stim.line = [];
stim.time = [];
stim.bisect = [];
response = struct;
time = struct;
data = struct;

% Dimensions in cm
sd = 57; %distance from laptop
sx = 48; %size of screen
sy = 27;

% Size of stimuli
stim.line.width = 0.3; %line width in degrees
stimlinewdth_cm = tan(pi*stim.line.width/360)*2*sd; %in cm
stim.line.offset = 0; %offset left or right by 5 degrees
stimlineoff_cm = tan(pi*stim.line.offset/360)*2*sd;
% Different linelengths
stim.line.lengths = [10, 20, 30];
stimlinelength_cm = tan(pi*stim.line.lengths/360)*2*sd; 
stim.line.dist = 0.75; %distance of line away from centre (on above, one below)
stimlinedist_cm = tan(pi*stim.line.dist/360)*2*sd; 

% Keyboard configuration
upkey = KbName('1'); 
downkey = KbName('2'); 

% Import data matrix
[Data,Text] = xlsread('TrialMatrix_LM'); %importing trial information for landmark tasks
stim.side = Data(:, ismember(Text, 'side')); %5deg left, 0, 5deg right
stim.size = Data(:, ismember(Text, 'size'));
stim.lshift = Data(:, ismember(Text, 'left shift'));
stim.rshift = Data(:, ismember(Text, 'right shift'));

% Deciding which stimuli to place  shift
% if 1 = shift on top line bottom remains equal, if 2 = shift on bottom
% line top remains equal
upper = zeros(135,1); upper = upper+1; lower = upper+1;
position = [upper; lower];
stim.pos = position(randperm(length(position)));

% Stimulus shift converting to pixels
% Left shift
stim.lshift = stim.lshift/10;
stim.rshift = stim.rshift/10;
for i = 1:length(stim.lshift)
    lshiftcm(i,:) = tan(pi*stim.lshift(i)/360)*2*sd;
    rshiftcm(i,:) = tan(pi*stim.rshift(i)/360)*2*sd;
end

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

%
for i = 1:nrtrials
    if stim.pos(i) == 1
       response.LHS(i) = upkey;
       response.RHS(i) = downkey;
    elseif stim.pos(i) == 2
       response.LHS(i) = downkey;
       response.RHS(i) = upkey;
    else 
       response.corr(i) = 5; %random number, can't be 0 because 0=no response, so would be corr if no response
    end
end 

        
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
HideCursor;

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
stim.line.lengthspix = stimlinelength_cm*DisplayXSize/sx; %pixels, 
stim.line.widthpix = stimlinewdth_cm*DisplayYSize/sy; 
%stim.dim.bisectlength.pix = stim.dim.bisectlength.cm*DisplayYSize/sy; %verticle so along y axis
stim.line.offsetpix = stimlineoff_cm*DisplayXSize/sx; 
stim.line.distpix = stimlinedist_cm*DisplayYSize/sy; 
%stim.dim.dist.pix = stim.dim.dist.cm*DisplayXSize/sx;

for i = 1:length(stim.lshift)
    stim.lshiftpix(i,:) = lshiftcm(i)*DisplayXSize/sx; 
    stim.rshiftpix(i,:) = rshiftcm(i)*DisplayXSize/sx; 
end

%% Starting experiment
% Welcome screen
Screen('FillRect', window, grey);
Screen('TextStyle', window, 0); %normal
Screen('TextFont', window, 'Lucida Console');
text = 'Press a button to start landmarks task';
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
    % Coordinates for horizontal line
    liney1 = midY-stim.line.distpix; %upperline
    liney2 = midY+stim.line.distpix; %lowerline
    linex1 = midX-linelength/2;
    linex2 = midX+linelength/2;
    
    % Stimulus side (from midline)
    switch stim.side(i)
        case 0
            offset = 0; %also the location for the middle of the line
        case 1
            offset = -stim.line.offsetpix; %offset left, subtract from x
        case 2
            offset = stim.line.offsetpix; %offset right, subtract from x
    end
    
    %% Present stim
    text = '+';
    width = RectWidth(Screen ('TextBounds', window, text));
    height = RectHeight(Screen ('TextBounds', window, text));
    Screen('Flip', window);
    Screen('DrawText', window, text, midX+offset, midY-height/2, black); %fixation cross at center
    Screen('Flip', window);
    WaitSecs(stim.time.fix)

        % Draw two lines
    if stim.pos(i) == 1
        % Top line
        Screen('DrawLine', window, black, linex1+offset-stim.lshiftpix(i), liney1-(stim.line.widthpix/2),...
            midX+offset, liney1-(stim.line.widthpix/2), (stim.line.widthpix/2)); %line 1 - top left
        Screen('DrawLine', window, white ,midX+offset, liney1-(stim.line.widthpix/2),...
            linex2+offset-stim.lshiftpix(i), liney1-(stim.line.widthpix/2), (stim.line.widthpix/2)); % line 2 - top right
        Screen('DrawLine', window, white, linex1+offset-stim.lshiftpix(i), liney1,...
            midX+offset, liney1, (stim.line.widthpix/2)); %line 3 - bottom left
        Screen('DrawLine', window, black ,midX+offset, liney1,...
            linex2+offset-stim.lshiftpix(i), liney1, (stim.line.widthpix/2)); % line 2 - bottom right
        % Bottom line
        Screen('DrawLine', window, white, linex1+offset+stim.rshiftpix(i), liney2-(stim.line.widthpix/2),...
            midX+offset, liney2-(stim.line.widthpix/2), (stim.line.widthpix/2)); %line 1 - top left
        Screen('DrawLine', window, black ,midX+offset, liney2-(stim.line.widthpix/2),...
            linex2+offset+stim.rshiftpix(i), liney2-(stim.line.widthpix/2), (stim.line.widthpix/2)); % line 2 - top right
        Screen('DrawLine', window, black, linex1+offset+stim.rshiftpix(i), liney2,...
            midX+offset, liney2, (stim.line.widthpix/2)); %line 3 - bottom left
        Screen('DrawLine', window, white ,midX+offset, liney2,...
            linex2+offset+stim.rshiftpix(i), liney2, (stim.line.widthpix/2)); % line 2 - bottom right
    elseif stim.pos(i) == 2        
        % Top line
        Screen('DrawLine', window, black, linex1+offset+stim.rshiftpix(i), liney1-(stim.line.widthpix/2),...
            midX+offset, liney1-(stim.line.widthpix/2), (stim.line.widthpix/2)); %line 1 - top left
        Screen('DrawLine', window, white, midX+offset, liney1-(stim.line.widthpix/2),...
            linex2+offset+stim.rshiftpix(i), liney1-(stim.line.widthpix/2), (stim.line.widthpix/2)); % line 2 - top right
        Screen('DrawLine', window, white, linex1+offset+stim.rshiftpix(i), liney1,...
            midX+offset, liney1, (stim.line.widthpix/2)); %line 3 - bottom left
        Screen('DrawLine', window, black ,midX+offset, liney1,...
            linex2+offset+stim.rshiftpix(i), liney1, (stim.line.widthpix/2)); % line 2 - bottom right
        % Bottom line
        Screen('DrawLine', window, white, linex1+offset-stim.lshiftpix(i), liney2-(stim.line.widthpix/2),...
            midX+offset, liney2-(stim.line.widthpix/2), (stim.line.widthpix/2)); %line 1 - top left
        Screen('DrawLine', window, black, midX+offset, liney2-(stim.line.widthpix/2),...
            linex2+offset-stim.lshiftpix(i), liney2-(stim.line.widthpix/2), (stim.line.widthpix/2)); % line 2 - top right
        Screen('DrawLine', window, black, linex1+offset-stim.lshiftpix(i), liney2,...
            midX+offset, liney2, (stim.line.widthpix/2)); %line 3 - bottom left
        Screen('DrawLine', window, white ,midX+offset, liney2,...
            linex2+offset-stim.lshiftpix(i), liney2, (stim.line.widthpix/2)); % line 2 - bottom right
    end
    
    Screen('Flip', window);
    time.target.start(i) = GetSecs;
    WaitSecs(stim.time.line)
    Screen('Flip',window);
    time0=GetSecs; %for response
    %KbStrokeWait;
    
    
    % Wait for a response
    flag2=0;
    while (GetSecs-time0 < stim.time.response )&& flag2==0
        [keyIsDown, secs, keyCode] = KbCheck; %check keyboard
        if keyIsDown %if any keypress
            if any(ismember(KbName(keyCode),[KbName(upkey), KbName(downkey)]))
                flag2 = 1; %
                time.keypress(i) = GetSecs;     % time for the first key press
                response.RT(i) = time.keypress(i) - time.target.start(i); %RT for the first key press
                key = KbName(keyCode);   % add the key name for the first key press

                Screen('Flip', window);
                time.target.end(i) = GetSecs;
            end
        end
    end
    
    if flag2==1 %some response
        if any(ismember(KbName(keyCode),KbName(downkey)))
            response.key(i) = downkey;%this is the code of the key pressed
        end

        if any(ismember(KbName(keyCode),KbName(upkey)))
            response.key(i) = upkey;%this is the code of the key pressed
        end 
    end
        
    if flag2==0 %no response
       response.key(i)=0;
       RT(i)=NaN;
       timekeypress(i)=NaN;
       WaitSecs(stim.time.response);  
       Screen('Flip', window);
       time.target.end(i) = GetSecs;
    end

end
time.experiment.end = GetSecs;

% Thank you screen
    thankYouText = 'End of landmarks, thank you!';
    height=RectHeight(Screen ('TextBounds', window, thankYouText));
    width=RectWidth(Screen ('TextBounds', window, thankYouText));
    Screen('DrawText', window, thankYouText, DisplayXSize/2 - width/2, DisplayYSize/2, white);
    Screen('Flip', window);
    WaitSecs (2);
    Screen('Flip', window);
    Screen('CloseAll');

%% Analysing
% Results
% Getting expected response
for i = 1:length(response.key)
    response.LHSlong(i) = isequal(response.LHS(i), response.key(i));
    response.RHSlong(i) = isequal(response.RHS(i), response.key(i));
end
data.RT = response.RT*1000; %transforming to ms

time.target.dur = time.target.end - time.target.start;
time.experiment.dur = time.experiment.end - time.experiment.start;
response.means.LHSlong = mean(response.LHSlong)*100;
response.means.RHSlong = mean(response.RHSlong)*100;

%% Saving data
catch
    rethrow(lasterror) 
    sca
    save(matfilename, 'stim', 'response', 'time', 'data')
end
sca
save(matfilename, 'stim', 'response', 'time', 'data')