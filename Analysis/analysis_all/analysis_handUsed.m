%% AG.Mitchell 30.04.20
% Script to investigate whether hand used to complete MLB and TRB task
% affects bisection error
% And whether reliability across time is increased if error is from one
% hand only 
% Addressing reviewer 2 comments in Neuropsychologia re-submission (1.1)

%% Setting paths
nParticipants = [1:19,21:24,26:30];
for p = 1:length(nParticipants)
    % Directory
    dirBias = ('M:\Alex_Files\Experiments\Bias'); %subject to change depending on where you analyse
    
    ppID = sprintf('P%0*d',2,nParticipants(p)); %for use when navigating files
    visualFilename = sprintf('%s_visualanalysisStart.mat', ppID);
    tactileFilename = sprintf('%s_tactileanalysis.mat', ppID);
    nSessions = 1:4; %vector number of sessions each participant does
    
    % individual participant directory, analysis
    dirPP = [dirBias filesep 'Data' filesep ppID filesep 'Analysis']; %participant directory
    dirVis = [dirPP filesep 'Visual'];
    dirTac = [dirPP filesep 'Tactile'];
    cd(dirPP)
    
%% Extracting data
    % for MLB - go to visual folder
    cd(dirVis)
    load(visualFilename)
    cd(dirTac)
    load(tactileFilename)
    
    %% Getting data per session, compiling into matrices
    %%%%%% reached here and keep going
    for i = 1:length(nSessions)
        mlbData.(sprintf('%s', session)).left_hand(:,p) = mlb.(sprintf('%s', session)).erro.left_hand(1);
    end
    
end
%% Group averages
%% Plotting data

%% Group averages
%% Plotting data

%% Alpha
%% Organising data for ANOVA