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
mkdir Tactile
dirVis = [dirAna 'Tactile' filesep];

trb = struct;
tr2 = struct;

%% Importing data for both tasks
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));
    dirSess = [dirPP filesep session filesep]; %current session pathway
    cd(dirSess); %directing to current session folder

    %% TRB task
    [Data,Text] = xlsread('TrialMatrix_TRB'); %importing trial information for landmark tasks
    size = Data(:, ismember(Text(1,:), 'size'));
    side = Data(:, ismember(Text(1,:), 'shift')); %2deg left = 1, middle = 0, 2deg right = 2
    response = Data(:, ismember(Text(1,:), 'response')); %in mm, relative to length of line
    % Finding error 
    size_mm = size*100; %finding actual line length in mm to calculate error
    error = response - size_mm/2; %calculating data error in mm
    trb.(sprintf('%s', session)).data.error = error;
    % data matrix for later analyses
    trb.(sprintf('%s', session)).data.matrix = [size, side, response, error];
    trb.(sprintf('%s', session)).error(1) = nanmean(error); %average results in mm
    trb.(sprintf('%s', session)).error(2) = nanstd(error);
    
    %% TR2 task
    [Data,Text] = xlsread('TrialMatrix_TR2'); %importing trial information for landmark tasks
    size = Data(:, ismember(Text(1,:), 'size'));
    shift = Data(:, ismember(Text(1,:), 'shift')); %shift in mm (0,2,4,6,8,1(=10))
    line = Data(:, ismember(Text(1,:), 'line')); %indicating the line that should be shifted left (the other right, 1 = top, 2 = bottom)
    response = Data(:, ismember(Text(1,:), 'response')); %top line = 1, bottom line = 2
    % Making matrix
    tr2.(sprintf('%s', session)).data.matrix = [size, shift, line, response];
end

%% TRB analysis
% Grouping into line length
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));
    cd(dirSess); %directing to current session folder
    % 10 cm line
    trb.(sprintf('%s', session)).line1mat = ...
        trb.(sprintf('%s', session)).matrix(find(trb.(sprintf('%s', session)).matrix(:,1)== 1),:);
    % 20 cm line
    trb.(sprintf('%s', session)).line2mat = ...
        trb.(sprintf('%s', session)).matrix(find(trb.(sprintf('%s', session)).matrix(:,1)== 2),:);
    % 30cm line
    trb.(sprintf('%s', session)).line3mat = ...
        trb.(sprintf('%s', session)).matrix(find(trb.(sprintf('%s', session)).matrix(:,1)== 3),:);
end
%% TRB plots
%% TR2 analysis
%% TR2 plots
%% Close and save

