%% Analysis for Tactile tasks - AG.Mitchell 01.10.18
% Loads in the data for the task and averages the response across all
% variables
% Then calculates bias
% Then calculates bias relative to line length and shift into L/R hemifield
clear all

%% File paths
ppID = input('Participant ID? ', 's'); %for use when navigating files
matfilename = sprintf('%s_tactileanalysis.mat', ppID);
nSessions = 1:4; %vector number of sessions each participant does
% Directory
dirBias = ('C:\Users\Experimenter\Documents\Experiments2018\Bias'); %subject to change depending on where you analyse
dirPP = [dirBias filesep ppID]; %participant directory
% Making new anaysis folder for saving
cd(dirPP)
mkdir Analysis;
dirAna = [dirPP filesep 'Analysis' filesep];
cd(dirAna)
mkdir Visual
dirVis = [dirAna 'Visual' filesep];

trb = struct;
tr2 = struct;

%% Importing data for both tasks
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));
    dirSess = [dirPP filesep session filesep]; %current session pathway
    cd(dirSess); %directing to current session folder

    %% TRB task
    [Data,Text] = xlsread('TrialMatrix_TRB'); %importing trial information for landmark tasks
    trb.(sprintf('%s', session)).data.size = Data(:, ismember(Text(1,:), 'size'));
    trb.(sprintf('%s', session)).data.side = Data(:, ismember(Text(1,:), 'shift')); %2deg left = 1, middle = 0, 2deg right = 2
    trb.(sprintf('%s', session)).data.response = Data(:, ismember(Text(1,:), 'response')); %in mm, relative to length of line
    % Finding error 
    trb.(sprintf('%s', session)).data.size = data.size*100; %finding actual line length in mm to calculate error
    trb.(sprintf('%s', session)).data.error = data.response - data.size/2; %calculating data error in mm
    % data matrix for later analyses
    trb.(sprintf('%s', session)).data.matrix = [data.size, data.side, data.response, data.error];
    trb.(sprintf('%s', session)).error(1) = nanmean(data.error); %average results in mm
    trb.(sprintf('%s', session)).error(2) = nanstd(data.error);
    
    %% TR2 task
    [Data,Text] = xlsread('TrialMatrix_TR2'); %importing trial information for landmark tasks
    tr2.(sprintf('%s', session)).data.size = Data(:, ismember(Text(1,:), 'size'));
    tr2.(sprintf('%s', session)).data.shift = Data(:, ismember(Text(1,:), 'shift')); %shift in mm (0,2,4,6,8,1(=10))
    tr2.(sprintf('%s', session)).data.line = Data(:, ismember(Text(1,:), 'line')); %indicating the line that should be shifted left (the other right, 1 = top, 2 = bottom)
    tr2.(sprintf('%s', session)).data.response = Data(:, ismember(Text(1,:), 'response')); %top line = 1, bottom line = 2
    % Making matrix
    tr2.(sprintf('%s', session)).data.matrix = [data.size, data.shift, data.line, data.response];
end


