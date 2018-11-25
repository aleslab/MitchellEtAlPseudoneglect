%% Bias study group-level analyis
%% AG.Mitchell - 24.11.18
% This script takes the data from the MLB, LM and TRB tasks for each
% participant
% Data gets collated in a manner that's translatable to most statistical
% platforms
% Takes proportion error from MLB and TRB, and the threshold to 50%
% (stimulus asymmetry when the participant thinks the line is bisected down
% the midline) for LM tasks and uses Cronbach's alpha to assess whether the
% data is consistent across sessions and modalities

%% Loading data for each session
nParticipants = 1:17;
allData = struct;
for p = 1:length(nParticipants)
    ppID = sprintf('P%0*d',2,nParticipants(p)); %for use when navigating files
    % Variables 
    nSessions = 1:4;
    visualFilename = sprintf('%s_visualanalysis.mat', ppID);
    tactileFilename = sprintf('%s_tactileanalysis.mat', ppID);
    % Directories
    % Directory
    dirBias = ('C:\Users\Experimenter\Documents\Experiments2018\Bias'); %subject to change depending on where you analyse
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
    
    %% Making structures and matrices for analysis
    for i = 1:length(nSessions)
        cd(dirPP)
        session = sprintf('Session%0*d',2,nSessions(i));
        
        % Saving data to matrix - one participant per row
        % LM task
        allData.(sprintf('%s', session)).lm.PSE(p,:) = lm.pfits.(sprintf('%s', session)).stim50right; %PSE at 50% for each participant
        allData.(sprintf('%s', session)).lm.CIs(p,:) = lm.pfits.(sprintf('%s', session)).threshCI;
        % Calculating SD from confidence intervals
        CIs = lm.pfits.(sprintf('%s', session)).threshCI;
        lmSD = sqrt(length(nParticipants))*(CIs(2)-CIs(1))/3.92;
        allData.(sprintf('%s', session)).lm.SDs(p,:) = lmSD;
        
        % MLB task
        allData.(sprintf('%s', session)).mlb.error(p,:) = mlb.(sprintf('%s', session)).error.mean(1);
        allData.(sprintf('%s', session)).mlb.errorstd(p,:) = mlb.(sprintf('%s', session)).error.mean(2);
        % Averaging line proportion error (cause I didn't do this
        % previously)...
        propError = [mlb.(sprintf('%s', session)).proportionError.line1(1), mlb.(sprintf('%s', session)).proportionError.line2(1),...
            mlb.(sprintf('%s', session)).proportionError.line3(1)];
        mlb.(sprintf('%s', session)).proportionError.mean(1) = mean(propError);
        mlb.(sprintf('%s', session)).proportionError.mean(2) = std(propError);
        % Saving to the alldata structure
        allData.(sprintf('%s', session)).mlb.proportionError(p,:) = mean(propError);
        allData.(sprintf('%s', session)).mlb.proportionErrorStd(p,:) = std(propError);
        
        % Same for TRB
        error = [trb.(sprintf('%s', session)).error.line1(1), trb.(sprintf('%s', session)).error.line2(1),...
            trb.(sprintf('%s', session)).error.line3(1)];
        trb.(sprintf('%s', session)).error.mean(1) = mean(error);
        trb.(sprintf('%s', session)).error.mean(2) = std(error);
        % Saving to all data structure
        allData.(sprintf('%s', session)).trb.error(p,:) = trb.(sprintf('%s', session)).error.mean(1);
        allData.(sprintf('%s', session)).trb.errorstd(p,:) = trb.(sprintf('%s', session)).error.mean(2);
        % Averaging line proportion error (cause I didn't do this
        % previously)...
        propError = [trb.(sprintf('%s', session)).proportionError.line1(1), trb.(sprintf('%s', session)).proportionError.line2(1),...
            trb.(sprintf('%s', session)).proportionError.line3(1)];
        trb.(sprintf('%s', session)).proportionError.mean(1) = mean(propError);
        trb.(sprintf('%s', session)).proportionError.mean(2) = std(propError);
        % Saving to the alldata structure
        allData.(sprintf('%s', session)).trb.proportionError(p,:) = mean(propError);
        allData.(sprintf('%s', session)).trb.proportionErrorStd(p,:) = std(propError);
        
        %% Saving data for  all sessions
        % Landmarks
        sessionlm =  allData.(sprintf('%s', session)).lm.PSE(p,1);
        sessionlmCI = allData.(sprintf('%s', session)).lm.CIs(p,:);
        sessionlmSD = allData.(sprintf('%s', session)).lm.SDs(p,1);
        allData.sessions.lmPSE(p,i) = sessionlm;
        allData.sessions.CIlow(p,i) = sessionlmCI(1);
        allData.sessions.CIhigh(p,i) = sessionlmCI(2);
        allData.sessions.lmStd(p,i) = sessionlmSD;
        % Manual line bisection
        sessionmlb = allData.(sprintf('%s', session)).mlb.proportionError(p,1);
        sessionmlbSD = allData.(sprintf('%s', session)).mlb.proportionErrorStd(p,1);
        allData.sessions.mlbProp(p,i) = sessionmlb;
        allData.sessions.mlbStd(p,i) = sessionmlbSD;
        % Tactile rod bisection
        sessiontrb = allData.(sprintf('%s', session)).trb.proportionError(p,1);
        sessiontrbSD = allData.(sprintf('%s', session)).trb.proportionErrorStd(p,1);
        allData.sessions.trbProp(p,i) = sessiontrb;
        allData.sessions.trbStd(p,i) = sessiontrbSD;
    end
    
    %% Saving data for all modalities
    % Modalities matrix: lm, mlb, trb; sd matrix the same
    % SD for LM calculated through confidence intervals
    % Landmarks
    modlm = mean(allData.sessions.lmPSE(p,:));
    modlmSD = std(allData.sessions.lmPSE(p,:));
    allData.modalities.data(p,1) = modlm; allData.modalities.sds(p,1) = modlmSD;
    % Manual line bisection
    modmlb = mean(allData.sessions.mlbProp(p,:));
    modmlbSD = std(allData.sessions.mlbProp(p,:));
    allData.modalities.data(p,2) = modmlb; allData.modalities.sds(p,2) = modmlbSD;
    % Tactile rod bisection
    modtrb = mean(allData.sessions.trbProp(p,:));
    modtrbSD = std(allData.sessions.trbProp(p,:));
    allData.modalities.data(p,3) = modtrb; allData.modalities.sds(p,3) = modtrbSD;
end

%% Next step: removing outliers from the analysis
% P11 I'm looking at you!