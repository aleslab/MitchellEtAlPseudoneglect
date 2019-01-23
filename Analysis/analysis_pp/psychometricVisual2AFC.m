%% AG. Mitchell - 23.01.19
%% Analysis for landmarks 2AFC task
% Takes the data from the landmarks 2AFC task and calculates percentage of
% 'top line longer' responses for each left/right shift.
% Calculates psychometric curve (cumulative normal) based on this data

%% Load data
%% Edited demo to fit bias study data (LM task) AG.Mitchell 22.11.18
clear all;      %Clear all existing variables from memory

tic
% Load in participant data
%nParticipants = [1:19,21,22,24];
nParticipants = 2; %for testing
for p = 1:length(nParticipants)  
    ppID = sprintf('P%0*d',2,nParticipants(p));
    %ppID = input('Participant ID? ', 's'); %for use when navigating files
    visfilename = sprintf('%s_visualanalysisStart.mat', ppID);
    matfilename = sprintf('%s_visualanalysis2AFC.mat', ppID);
    nSessions = 1:4; %vector number of sessions each participant does
    % Directory
    dirBias = ('C:\Users\Experimenter\Documents\Experiments2018\Bias'); %subject to change depending on where you analyse
    dirPP = [dirBias filesep ppID]; %participant directory
    % Navigating to analysis folder
    cd(dirPP)
    dirAna = [dirPP filesep 'Analysis' filesep];
    cd(dirAna)
    dirVis = [dirAna 'Visual' filesep];
    cd(dirVis)
    load(visfilename)

%% Calculate proportion response
%% Psychometric curve fitting
end