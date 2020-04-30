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
    matfilename = sprintf('%s_visualanalysisStart.mat', ppID);
    nSessions = 1:4; %vector number of sessions each participant does
    
    % individual participant directory, analysis
    dirPP = [dirBias filesep 'Data' filesep ppID filesep 'Analysis']; %participant directory
    dirVis = [dirPP filesep 'Visual'];
    dirTac = [dirPP filesep 'Tactile'];
    cd(dirPP)
    
%% MLB Extracting data
    % for MLB - go to visual folder
    
%% TRB Extracting data
end
%% Group averages
%% Plotting data

%% Group averages
%% Plotting data

%% Alpha
%% Organising data for ANOVA