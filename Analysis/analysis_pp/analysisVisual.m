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
    mlb.data.(sprintf('%s', session)) = data;
    
    % Landmarks
    lmName = dir([dirSess 'LM_*.mat']); %getting file details for MLB data
    load(lmName.name);    
    lm.data.(sprintf('%s', session)) = data;
    
    % Landmarks 2AFC
    lm2Name= dir([dirSess 'LM2afc_*.mat']); %getting file details for MLB data
    load(lm2Name.name);    
    lm2.data.(sprintf('%s', session)) = data;
end

%% Analyse MLB data
% Grouping into line length
for i = 1:length(nSessions)
end

%% Analyse LM data
%% Analyse LM data
%% Close and save