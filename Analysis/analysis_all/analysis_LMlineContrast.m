%% AG.Mitchell - 27.05.20
%% Analysis script looking at differences in response to different line contrasts for the landmark task
% As requested by reviewer 1 in Neuropsychologia resubmission (1.1)
% Script runs through right side longer responses for each line contrast (
% top left = black; top left = white) to see if visual illusion may be
% biasing response 

%% Setting paths
clear all
nParticipants = [1:19,21:24,26:30];
% Directory
dirBias = ('M:\Alex_Files\Experiments\Bias'); %subject to change depending on where you analyse
dirAna = ('M:\Alex_Files\Experiments\Bias\Analysis');

%% Looping through participants and getting data
% Average across sessions
% Individual sessions (1-4)
for p = 1:length(nParticipants)
    ppID = sprintf('P%0*d',2,nParticipants(p)); %for use when navigating files
    visualFilename = sprintf('%s_visualanalysisStart.mat', ppID);
    matfilename = ('LMcontrast_analysis.mat');
    
    nSessions = 1:4; %vector number of sessions each participant does
    
    % individual participant directory
    dirPP = [dirBias filesep 'Data' filesep ppID filesep 'Analysis']; %participant directory
    dirVis = [dirPP filesep 'Visual'];
    cd(dirVis)
    % load file
    load(visualFilename)
    
    
end