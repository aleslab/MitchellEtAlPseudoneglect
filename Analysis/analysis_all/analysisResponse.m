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
%nParticipants = 3; %for testing
nSessions = [1:4];
allData = struct;
matfilename = ('ReliabilityAnalysis.mat');
for p = 1:length(nParticipants)
    ppID = sprintf('P%0*d',2,nParticipants(p)); %for use when navigating files
    % Variables 
    nSessions = 1:4;
    visualFilename = sprintf('%s_visualanalysis2AFC.mat', ppID);
    tactileFilename = sprintf('%s_tactileanalysis.mat', ppID);
    matFilename = ('responseAnalysis.mat');
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
    for i = 1:length(nSessions)
        cd(dirPP)
        session = sprintf('Session%0*d',2,nSessions(i));
        % loading data into session matrix each column is individual
        % participant
        allData.(sprintf('%s', session)).names{p} = ppID;
        % landmarks data
        allData.(sprintf('%s', session)).lm2.asymmetry = lm2psych.(sprintf('%s', session))(:,1);
        allData.(sprintf('%s', session)).lm2.allPP(:,(p)) = lm2psych.(sprintf('%s', session))(:,2);
        % tactile data
        allData.(sprintf('%s', session)).tr2.asymmetry = tr2psych.(sprintf('%s', session))(:,1);
        allData.(sprintf('%s', session)).tr2.allPP(:,(p)) = tr2psych.(sprintf('%s', session))(:,2);
    end
    
    % Matrices outside sessions- average across all sessions for all PP
    allData.(sprintf('%s', session)).names{p} = ppID; %adding participant ID to name
    % landmarks
    allData.allSessions.lm2.asymmetry = lm2psych.allSessions(:,1);
    allData.allSessions.lm2.allPP(:,(p)) = lm2psych.allSessions(:,2);
    % tactile
    allData.allSessions.tr2.asymmetry = tr2psych.allSessions(:,1);
    allData.allSessions.tr2.allPP(:,(p)) = tr2psych.allSessions(:,2);
end

%% Averaging across participants
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));
    % landmarks
    allData.(sprintf('%s', session)).lm2.av(:,1) = mean(allData.(sprintf('%s', session)).lm2.allPP, 2); %mean
    allData.(sprintf('%s', session)).lm2.av(:,2) = std(allData.(sprintf('%s', session)).lm2.allPP, 0,2); 
    % tactile
    allData.(sprintf('%s', session)).tr2.av(:,1) = mean(allData.(sprintf('%s', session)).tr2.allPP, 2); %mean
    allData.(sprintf('%s', session)).tr2.av(:,2) = std(allData.(sprintf('%s', session)).tr2.allPP, 0,2); 
end
% across all sessions
allData.allSessions.lm2.av(:,1) = mean(allData.allSessions.lm2.allPP, 2); %mean 
allData.allSessions.lm2.av(:,2) = std(allData.allSessions.lm2.allPP, 0,2); %std

allData.allSessions.tr2.av(:,1) = mean(allData.allSessions.tr2.allPP, 2); %mean 
allData.allSessions.tr2.av(:,2) = std(allData.allSessions.tr2.allPP, 0,2); %std

%% Plotting mean data
% Make a plot that runs from -10 to 0 asymmetry as both lines shifted by
% the same amount
pdfFileName = ('responseAveragePP_sessions.mat');
figure()
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));
    % plotting
    M = plot(allData.allSessions.lm2.asymmetry(1:6), allData.(sprintf('%s', session)).lm2.av(1:6,1), ...
        'LineWidth', 1.5);
    hold on
end

%% Making Dakin plot
%% Binomial tests 
%% Plot Binomial