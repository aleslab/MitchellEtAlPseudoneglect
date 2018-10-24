%% A.G. Mitchell 22.10.18
%% Imports MLB data from each session for analysis
% This is the analysis script run on each participant
% Takes the average error from each session and calculates total
% participant bias
% Also important to maintain mean biases in each session for comparison in
% the second (across time) hypothesis
% Places data in a format readable for the overall analysis

clear all

%% Variables
ppID = input('Participant ID? ', 's'); %for use when navigating files
matfilename = sprintf('MLB_analysis%s.mat', ppID);
nSessions = 1:3; %vector number of sessions each participant does
% Directory
dirBias = ('C:\Users\Experimenter\Documents\Experiments2018\Bias'); %subject to change depending on where you analyse
dirPP = [dirBias filesep ppID filesep]; %participant directory

% Load data for all trials
cd(dirPP)
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));
    dirSess = [dirPP filesep session filesep]; %current session pathway
    cd(dirSess); %directing to current session folder
    
    % Manual line bisection
    mlbName = dir([dirSess 'MLB_*.mat']); %getting file details for MLB data
    load(mlbName.name);    
    mlb.(sprintf('%s', session)) = data;
    
    % Landmarks
    lmName = dir([dirSess 'LM_*.mat']); %getting file details for MLB data
    load(lmName.name);    
    lm.(sprintf('%s', session)) = data;
    
    % Landmarks 2AFC
    lm2Name= dir([dirSess 'LM2afc_*.mat']); %getting file details for MLB data
    load(lm2Name.name);    
    lm2.(sprintf('%s', session)) = data;
end

%% Analyse MLB data
% Grouping into line length
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));
    % 10 cm line
    mlb.(sprintf('%s', session)).line1mat = ...
        mlb.(sprintf('%s', session)).matrix(find(mlb.(sprintf('%s', session)).matrix(:,1)== 1),:);
    % 20 cm line
    mlb.(sprintf('%s', session)).line2mat = ...
        mlb.(sprintf('%s', session)).matrix(find(mlb.(sprintf('%s', session)).matrix(:,1)== 2),:);
    % 30cm line
    mlb.(sprintf('%s', session)).line3mat = ...
        mlb.(sprintf('%s', session)).matrix(find(mlb.(sprintf('%s', session)).matrix(:,1)== 3),:);
    
    % Average and std line length for each session
    avL1 = nanmean(mlb.(sprintf('%s', session)).line1mat(:,6)); %10 cm line
    stdL1 = nanstd(mlb.(sprintf('%s', session)).line1mat(:,6));
    mlb.(sprintf('%s', session)).error.line1 = [avL1, stdL1];
    avL2 = nanmean(mlb.(sprintf('%s', session)).line2mat(:,6)); %20 cm line
    stdL2 = nanstd(mlb.(sprintf('%s', session)).line2mat(:,6));
    mlb.(sprintf('%s', session)).error.line2 = [avL2, stdL2];
    avL3 = nanmean(mlb.(sprintf('%s', session)).line3mat(:,6)); %30 cm line
    stdL3 = nanstd(mlb.(sprintf('%s', session)).line3mat(:,6));
    mlb.(sprintf('%s', session)).error.line3 = [avL3, stdL3];
end



%% Analyse LM data 
% Take percentage left-side longer responses for each shift in mm
% Grouping into line length
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));
    % 10 cm line
    lm.(sprintf('%s', session)).line1mat = ...
        lm.(sprintf('%s', session)).matrix(find(lm.(sprintf('%s', session)).matrix(:,1)== 1),:);
    % 20 cm line
    lm.(sprintf('%s', session)).line2mat = ...
        lm.(sprintf('%s', session)).matrix(find(lm.(sprintf('%s', session)).matrix(:,1)== 2),:);
    % 30cm line
    lm.(sprintf('%s', session)).line3mat = ...
        lm.(sprintf('%s', session)).matrix(find(lm.(sprintf('%s', session)).matrix(:,1)== 3),:);
end

%% Analyse LM data
% Take percentage top-line longer responses for each shift in mm (of the
% top line)
% Grouping into line length
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));
    % 10 cm line
    lm2.(sprintf('%s', session)).line1mat = ...
        lm2.(sprintf('%s', session)).matrix(find(lm2.(sprintf('%s', session)).matrix(:,1)== 1),:);
    % 20 cm line
    lm2.(sprintf('%s', session)).line2mat = ...
        lm2.(sprintf('%s', session)).matrix(find(lm2.(sprintf('%s', session)).matrix(:,1)== 2),:);
    % 30cm line
    lm2.(sprintf('%s', session)).line3mat = ...
        lm2.(sprintf('%s', session)).matrix(find(lm2.(sprintf('%s', session)).matrix(:,1)== 3),:);
end

%% Close and save