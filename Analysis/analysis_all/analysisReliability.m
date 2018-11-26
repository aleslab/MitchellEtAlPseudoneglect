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
matfilename = ('ReliabilityAnalysis.mat');
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

%% Removing outliers from the analysis
% Any outliers identified through MATLAB - more than three scaled mean
% absolute deviations away
% Landmarks
outlierLM = isoutlier(allData.sessions.lmPSE);
outlierSumLM = sum(outlierLM,2); %sums the number of outliers across all participants
% If the sum of outliers is >2, remove participant from the analysis
removeLM = find(outlierSumLM > 2); %identifies paticipants that need removing

% Manual line bisection
outlierMLB = isoutlier(allData.sessions.mlbProp);
outlierSumMLB = sum(outlierMLB,2); %sums the number of outliers across all participants
% If the sum of outliers is >2, remove participant from the analysis
removeMLB = find(outlierSumMLB > 2); %identifies paticipants that need removing

% Tactile rod bisection
outlierTRB = isoutlier(allData.sessions.trbProp);
outlierSumTRB = sum(outlierTRB,2); %sums the number of outliers across all participants
% If the sum of outliers is >2, remove participant from the analysis
removeTRB = find(outlierSumTRB > 2); %identifies paticipants that need removing

removeAll = [removeLM; removeMLB; removeTRB]; 
removeAll = sort(removeAll); %sorting so in participant order

% Removing the participants found to be outliers
% Modalities and sessions
for r = 1:length(removeAll)
    % Identifying the outlier
    allData.sessions.lmPSE(removeAll(r),:) = NaN;
    allData.sessions.lmCIlow(removeAll(r),:) = NaN;
    allData.sessions.lmCIhigh(removeAll(r),:) = NaN;
    allData.sessions.lmStd(removeAll(r),:) = NaN;
    allData.sessions.mlbProp(removeAll(r),:) = NaN;
    allData.sessions.mlbPropStd(removeAll(r),:) = NaN;
    allData.sessions.trbProp(removeAll(r),:) = NaN;
    allData.sessions.trbPropStd(removeAll(r),:) = NaN;
    
    % Modalities across sessions
    allData.modalities.data(removeAll(r),:) = NaN;
    allData.modalities.sds(removeAll(r),:) = NaN;
end
% Removing outlier from matrix - modalities and sessions
allData.sessions.lmPSE = allData.sessions.lmPSE(find(allData.sessions.lmPSE(:,1) > 0.00001 | ...
    allData.sessions.lmPSE(:,1) < 0.00001),:);
allData.sessions.lmCIhigh= allData.sessions.lmCIhigh(find(allData.sessions.lmCIhigh(:,1) > 0.00001 | ...
    allData.sessions.lmCIhigh(:,1) < 0.00001),:);
allData.sessions.lmCIlow= allData.sessions.lmCIlow(find(allData.sessions.lmCIlow(:,1) > 0.00001 | ...
    allData.sessions.lmCIlow(:,1) < 0.00001),:);
allData.sessions.lmStd= allData.sessions.lmStd(find(allData.sessions.lmStd(:,1) > 0.00001 | ...
    allData.sessions.lmStd(:,1) < 0.00001),:);
allData.sessions.mlbProp= allData.sessions.mlbProp(find(allData.sessions.mlbProp(:,1) > 0.00001 | ...
    allData.sessions.mlbProp(:,1) < 0.00001),:);
allData.sessions.mlbPropStd= allData.sessions.mlbPropStd(find(allData.sessions.mlbPropStd(:,1) > 0.00001 | ...
    allData.sessions.mlbPropStd(:,1) < 0.00001),:);
allData.sessions.trbProp= allData.sessions.trbProp(find(allData.sessions.trbProp(:,1) > 0.00001 | ...
    allData.sessions.trbProp(:,1) < 0.00001),:);
allData.sessions.trbPropStd= allData.sessions.trbPropStd(find(allData.sessions.trbPropStd(:,1) > 0.00001 | ...
    allData.sessions.trbPropStd(:,1) < 0.00001),:);

% Modalities across sessions
allData.modalities.data = allData.modalities.data(find(allData.modalities.data(:,1) > 0.00001 | ...
    allData.modalities.data(:,1) < 0.00001),:);
allData.modalities.sds = allData.modalities.sds(find(allData.modalities.sds(:,1) > 0.00001 | ...
    allData.modalities.sds(:,1) < 0.00001),:);

%% Calculating group means
% General rule of thumb: 1st column = mean, 2nd column - sd
% Across sessions
% Mean for each
allData.sessions.means.lmPSE = mean(allData.sessions.lmPSE);
allData.sessions.means.mlbProp = mean(allData.sessions.mlbProp);
allData.sessions.means.trbProp = mean(allData.sessions.trbProp);
% Std for each
allData.sessions.sds.lmPSE = std(allData.sessions.lmPSE);
allData.sessions.sds.mlbProp = std(allData.sessions.mlbProp);
allData.sessions.sds.trbProp = std(allData.sessions.trbProp);

% Total for each modality
allData.means.lmPSE(1) = mean(allData.sessions.means.lmPSE);
allData.means.lmPSE(2) = std(allData.sessions.means.lmPSE);
allData.means.mlbProp(1) = mean(allData.sessions.means.mlbProp);
allData.means.mlbProp(2) = std(allData.sessions.means.mlbProp);
allData.means.trbProp(1) = mean(allData.sessions.means.trbProp);
allData.means.trbProp(2) = std(allData.sessions.means.trbProp);

%Total
allData.means.tot(1) = mean([allData.means.lmPSE(1), allData.means.mlbProp(1), ...
    allData.means.trbProp(1)]);
allData.means.tot(2) = std([allData.means.lmPSE(1), allData.means.mlbProp(1), ...
    allData.means.trbProp(1)]);

%% Modalities plot
%% Sessions plot
%% Cronbach's alpha and other stats...

%% save and close
close all
cd(dirAnaAll)
save(matfilename, 'allData');
