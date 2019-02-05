%% AG.Mitchell - 05.02.19
%% All participant analysis for bias study results
% 'Response' hypothesis

% Getting average of all participant data and running binomial tests on
% this data
% Getting binomial distribution plots for each stimulus asymmetry
% Plotting the response version of the 'Dakin' plot summarising consistency
% across sessions for these tasks too (and overall participant bias)

clear all
%% Getting data
nParticipants = [1:19,21,22,24];
nSessions = [1:4];
allData = struct;
matfilename = ('ReliabilityAnalysis.mat');
for p = 1:length(nParticipants)
    ppID = sprintf('P%0*d',2,nParticipants(p)); %for use when navigating files
    % Variables 
    nSessions = 1:4;
    visualFilename = sprintf('%s_visualanalysis2AFC.mat.mat', ppID);
    tactileFilename = sprintf('%s_tactileanalysis.mat', ppID);
    % Directories
    % Directory
    dirBias = ('M:\Experiments\Bias'); %subject to change depending on where you analyse
    dirPP = [dirBias filesep ppID]; %participant directory
    dirAna = [dirPP filesep 'Analysis' filesep];
    dirVis = [dirAna 'Visual' filesep];
    dirTact = [dirAna 'Tactile' filesep];
    dirAnaAll = [dirBias filesep 'Analysis']; %directory for all analysis - here is where data should be saved from this file
    
    % Loading files for each participant visual and tactile data analysis
    cd(dirVis)
    load(visualFilename)
    cd(dirTact)
    load(tactileFilename)
    
    %% Make structures and matrices for analysis
    for j = 1:length(nSessions)
    end
    
    %% Plotting mean data
    %% Making Dakin plot
    %% Binomial tests 
    %% Plot Binomial
end