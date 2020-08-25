%% AG.Mitchell - 27.05.20
%% Analysis script looking at differences in response to different line contrasts for the landmark task
% As requested by reviewer 1 in Neuropsychologia resubmission (1.1)
% Script runs through right side longer responses for each line contrast (
% top left = black; top left = white) to see if visual illusion may be
% biasing response 

%% Setting paths
clear all
nParticipants = [1:10,12:16,18,21:24,26:30];
% Directory
% Getting directory
filePath = cd;
[dirBias, name, ext] = fileparts(filePath); %subject to change depending on where you analysis;
cd(dirBias)

%% Looping through participants and getting data
for p = 1:length(nParticipants)
    ppID = sprintf('P%0*d',2,nParticipants(p)); %for use when navigating files
    visualFilename = sprintf('%s_visualanalysisContrast.mat', ppID);
    matfilename = ('LMcontrast_analysis.mat');
    
    % individual participant directory
    dirPP = [dirBias filesep 'Data' filesep ppID filesep 'Analysis']; %participant directory
    dirVis = [dirPP filesep 'Visual'];
    cd(dirVis)
    % load file
    load(visualFilename)
    
    %% Extract data
    % Get PSE for each participant
    PSE_con1 = -(lm.pfits.con1.stim50right);
    PSE_con2 = -(lm.pfits.con2.stim50right);
    res.con1.PSE(p) = PSE_con1;
    res.con2.PSE(p) = PSE_con2;
end

%% Statistical analysis
% t-test comparing mean percentage right-side responses across all sessions
[res.all.h, res.all.p, res.all.ci, res.all.stats] = ttest(res.con1.PSE, res.con2.PSE);

%% Plotting
