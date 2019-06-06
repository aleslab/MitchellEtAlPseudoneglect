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
nParticipants = [1:19,21:24,26:30];
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
    
    %% Making structures and matrices for analysis
    for i = 1:length(nSessions)
        cd(dirPP)
        session = sprintf('Session%0*d',2,nSessions(i));
        
        % Saving data to matrix - one participant per row
        % LM task
        lmPSEdata = lm.pfits.(sprintf('%s', session)).stim50right; %PSE at 50% for each participant
        % Changing the sign of the landmarks PSE data to fit with more
        % negative more left' - as the more negative here shows that the
        % right side is viewed as longer. So now data will provide value of
        % how much longer does the participant view the biased side
        signChange = -(lmPSEdata);
        allData.(sprintf('%s', session)).lm.PSE(p,:) = signChange;
        allData.(sprintf('%s', session)).lm.CIs(p,:) = lm.pfits.(sprintf('%s', session)).threshCI; %adding CIs
        % Calculating SD from confidence intervals
        CIs = lm.pfits.(sprintf('%s', session)).threshCI;
        lmSD = sqrt(length(nParticipants))*(CIs(2)-CIs(1))/3.92;
        allData.(sprintf('%s', session)).lm.SDs(p,:) = lmSD;

        
        % MLB task
        allData.(sprintf('%s', session)).mlb.error(p,:) = (mlb.(sprintf('%s', session)).error.mean(1))*10; %converting back to mm
        allData.(sprintf('%s', session)).mlb.errorstd(p,:) = (mlb.(sprintf('%s', session)).error.mean(2))*10; %converting back to mm
        % Averaging line proportion error (cause I didn't do this
        % previously)...
        propError = [mlb.(sprintf('%s', session)).proportionError.line1(1), mlb.(sprintf('%s', session)).proportionError.line2(1),...
            mlb.(sprintf('%s', session)).proportionError.line3(1)];
        mlb.(sprintf('%s', session)).proportionError.mean(1) = nanmean(propError);
        mlb.(sprintf('%s', session)).proportionError.mean(2) = nanstd(propError);
        % Saving to the alldata structure
        allData.(sprintf('%s', session)).mlb.proportionError(p,:) = nanmean(propError);
        allData.(sprintf('%s', session)).mlb.proportionErrorStd(p,:) = nanstd(propError);
        
        % Same for TRB
        error = [trb.(sprintf('%s', session)).error.line1(1), trb.(sprintf('%s', session)).error.line2(1),...
            trb.(sprintf('%s', session)).error.line3(1)];
        trb.(sprintf('%s', session)).error.mean(1) = nanmean(error)*10; %converting back to mm;
        trb.(sprintf('%s', session)).error.mean(2) = nanstd(error); 
        % Saving to all data structure
        allData.(sprintf('%s', session)).trb.error(p,:) = trb.(sprintf('%s', session)).error.mean(1);
        allData.(sprintf('%s', session)).trb.errorstd(p,:) = trb.(sprintf('%s', session)).error.mean(2);
        % Averaging line proportion error (cause I didn't do this
        % previously)...
        propError = [trb.(sprintf('%s', session)).proportionError.line1(1), trb.(sprintf('%s', session)).proportionError.line2(1),...
            trb.(sprintf('%s', session)).proportionError.line3(1)];
        trb.(sprintf('%s', session)).proportionError.mean(1) = nanmean(propError);
        trb.(sprintf('%s', session)).proportionError.mean(2) = nanstd(propError);
        % Saving to the alldata structure
        allData.(sprintf('%s', session)).trb.proportionError(p,:) = nanmean(propError);
        allData.(sprintf('%s', session)).trb.proportionErrorStd(p,:) = nanstd(propError);
        
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
        sessionmlb = allData.(sprintf('%s', session)).mlb.error(p,1);
        sessionmlbSD = allData.(sprintf('%s', session)).mlb.error(p,1);
        allData.sessions.mlb(p,i) = sessionmlb;
        allData.sessions.mlbStd(p,i) = sessionmlbSD;
        % Proportion error
        sessionmlb = allData.(sprintf('%s', session)).mlb.proportionError(p,1);
        sessionmlbSD = allData.(sprintf('%s', session)).mlb.proportionErrorStd(p,1);
        allData.sessions.mlbProp(p,i) = sessionmlb;
        allData.sessions.mlbPropStd(p,i) = sessionmlbSD;
        
        % Tactile rod bisection
        sessiontrb = allData.(sprintf('%s', session)).trb.error(p,1);
        sessiontrbSD = allData.(sprintf('%s', session)).trb.error(p,1);
        allData.sessions.trb(p,i) = sessiontrb;
        allData.sessions.trbStd(p,i) = sessiontrbSD;
        % Proportion error
        sessiontrb = allData.(sprintf('%s', session)).trb.proportionError(p,1);
        sessiontrbSD = allData.(sprintf('%s', session)).trb.proportionErrorStd(p,1);
        allData.sessions.trbProp(p,i) = sessiontrb;
        allData.sessions.trbPropStd(p,i) = sessiontrbSD;
        
    end
    
    %% Saving data for all modalities
    % Modalities matrix: lm, mlb, trb; sd matrix the same
    % SD for LM calculated through confidence intervals
    % Landmarks
    modlm = nanmean(allData.sessions.lmPSE(p,:));
    modlmSD = nanstd(allData.sessions.lmPSE(p,:));
    allData.modalities.data(p,1) = modlm; allData.modalities.sds(p,1) = modlmSD;
    % Manual line bisection
    modmlb = nanmean(allData.sessions.mlb(p,:));
    modmlbSD = nanstd(allData.sessions.mlbProp(p,:));
    allData.modalities.data(p,2) = modmlb; allData.modalities.sds(p,2) = modmlbSD;
    % Tactile rod bisection
    modtrb = nanmean(allData.sessions.trb(p,:));
    modtrbSD = nanstd(allData.sessions.trb(p,:));
    allData.modalities.data(p,3) = modtrb; allData.modalities.sds(p,3) = modtrbSD;
    
    
    %% Mean for all PP
    % column 1 means, column 2 SDs
    allData.means.allPP(p,1) = nanmean([modlm, modmlb, modtrb]);
    allData.means.allPP(p,2) = nanstd([modlm, modmlb, modtrb]);
end

% Need sds across modalities for individual sessions
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));
    sessInfo = [allData.(sprintf('%s', session)).lm.PSE, allData.(sprintf('%s', session)).mlb.error, ...
        allData.(sprintf('%s', session)).trb.error];
    allData.sessions.(sprintf('%sall', session))(:,1) = nanmean(sessInfo,2);
    allData.sessions.(sprintf('%sall', session))(:,2) = nanstd(sessInfo,0,2);
end
%% Removing outliers from the analysis
% Any outliers identified through MATLAB - more than three scaled mean
% absolute deviations away
% Landmarks
outlierLM = isoutlier(allData.modalities.data(:,1));
% Manual line bisection
outlierMLB = isoutlier(allData.modalities.data(:,2));
% Tactile rod bisection
outlierTRB = isoutlier(allData.modalities.data(:,3));
outlierAll = [outlierLM, outlierMLB, outlierTRB]; outlierSum = sum(outlierAll,2);
% Participant 24 = curve fitted poorly and not identified through outlier removeal
% should be removed
outlierSum(24,1) = 1;

remove = find(outlierSum > 0.1); %identifies paticipants that need removing
remove = sort(remove); %sorting so in participant order

% Removing the participants found to be outliers
% Modalities and sessions
for r = 1:length(remove)
    % Identifying the outlier
    allData.sessions.lmPSE(remove(r),:) = NaN;
    allData.sessions.CIlow(remove(r),:) = NaN;
    allData.sessions.CIhigh(remove(r),:) = NaN;
    allData.sessions.lmStd(remove(r),:) = NaN;
    allData.sessions.mlb(remove(r),:) = NaN;
    allData.sessions.mlbStd(remove(r),:) = NaN;
    allData.sessions.trb(remove(r),:) = NaN;
    allData.sessions.trbStd(remove(r),:) = NaN;
    
    % Modalities within sessions
    for i = 1:length(nSessions)
        session = sprintf('Session%0*d',2,nSessions(i));
        allData.sessions.(sprintf('%sall', session))(remove(r),:) = NaN;
    end
    
    % Modalities across sessions
    allData.modalities.data(remove(r),:) = NaN;
    allData.modalities.sds(remove(r),:) = NaN;
    
    % All means
    allData.means.allPP(remove(r),:) = NaN;
end
% Removing outlier from matrix - modalities and sessions
allData.sessions.lmPSE = allData.sessions.lmPSE(isfinite(allData.sessions.lmPSE(:,1)),:);
allData.sessions.CIhigh= allData.sessions.CIhigh(isfinite(allData.sessions.CIhigh(:,1)),:);
allData.sessions.CIlow = allData.sessions.CIlow(isfinite(allData.sessions.CIlow(:,1)),:);
allData.sessions.lmStd= allData.sessions.lmStd(isfinite(allData.sessions.lmStd(:,1)),:);
allData.sessions.mlb= allData.sessions.mlb(isfinite(allData.sessions.mlb(:,1)),:);
allData.sessions.mlbStd= allData.sessions.mlbStd(isfinite(allData.sessions.mlbStd(:,1)),:);
allData.sessions.trb= allData.sessions.trb(isfinite(allData.sessions.trb(:,1)),:);
allData.sessions.trbStd= allData.sessions.trb(isfinite(allData.sessions.trb(:,1)),:);

% Modalities within sessions
allData.sessions.Session01all= allData.sessions.Session01all(isfinite(allData.sessions.Session01all(:,1)),:);
allData.sessions.Session02all= allData.sessions.Session02all(isfinite(allData.sessions.Session02all(:,1)),:);
allData.sessions.Session03all= allData.sessions.Session03all(isfinite(allData.sessions.Session03all(:,1)),:);
allData.sessions.Session04all= allData.sessions.Session04all(isfinite(allData.sessions.Session04all(:,1)),:);

% Modalities across sessions
allData.modalities.data = allData.modalities.data(isfinite(allData.modalities.data(:,1)),:);
allData.modalities.sds = allData.modalities.sds(isfinite(allData.modalities.sds(:,1)),:);

% All means
allData.means.allPP = allData.means.allPP(isfinite(allData.means.allPP(:,1)),:);

%% Calculating group means
% General rule of thumb: 1st column = mean, 2nd column - sd
% Across sessions
% Mean for each
allData.sessions.means.lmPSE = mean(allData.sessions.lmPSE);
allData.sessions.means.mlb = mean(allData.sessions.mlb);
allData.sessions.means.trb = mean(allData.sessions.trb);
% Std for each
allData.sessions.sds.lmPSE = std(allData.sessions.lmPSE);
allData.sessions.sds.mlb = std(allData.sessions.mlb);
allData.sessions.sds.trb = std(allData.sessions.trb);

% Total for each modality
allData.means.lmPSE(1) = mean(allData.sessions.means.lmPSE);
allData.means.lmPSE(2) = std(allData.sessions.means.lmPSE);
allData.means.mlb(1) = mean(allData.sessions.means.mlb);
allData.means.mlb(2) = std(allData.sessions.means.mlb);
allData.means.trb(1) = mean(allData.sessions.means.trb);
allData.means.trb(2) = std(allData.sessions.means.trb);

%Total
allData.means.tot(1) = mean([allData.means.lmPSE(1), allData.means.mlb(1), ...
    allData.means.trb(1)]);
allData.means.tot(2) = std([allData.means.lmPSE(1), allData.means.mlb(1), ...
    allData.means.trb(1)]);


%% Data plots
results = struct;
results.observers = 1:length(allData.means.allPP); %getting total number of observers

%% Modalities plot
% The idea for this plot was taken from the Dakin (2018, Vision Res) papert
% figure 2a
% Plot for mean bias in all participants
% First - need to order by bias (left -> right)
% Need to place all data in the same matrix to do this, then sort entire
% matrix by mean bias
results.plotting.modalities(:,1) = results.observers;
results.plotting.modalities(:,2) = allData.means.allPP(:,1); %mean all data
results.plotting.modalities(:,3) = allData.modalities.data(:,1); %lm
results.plotting.modalities(:,4) = allData.modalities.data(:,2); %mlb
results.plotting.modalities(:,5) = allData.modalities.data(:,3); %trb.
results.plotting.modalities(:,6) = allData.means.allPP(:,2); %std all data (across modalities)
% Getting confidence intervals for all tasks
CIall = (allData.means.allPP(:,2)/sqrt(length(results.observers)))*2.11; %s chosen here because that's the number of trials completed for each task in 1 session
results.plotting.modalities(:,7) = CIall;
% Sorting the matrix by mean bias
results.plotting.modalities = sortrows(results.plotting.modalities, 2);
results.plotting.modalities(:,8) = results.observers; %observers not sorted by mean bias for use with plotting - organisation of data

% standard deviation values for shading
SDpt5 = allData.means.tot(2)*0.5;
SD2 = allData.means.tot(2)*2;
SDall = results.plotting.modalities(:,6);
CIall = results.plotting.modalities(:,7);

cd(dirAnaAll);
% Making initial mean error (all modalities) plot
figure('units', 'centimeters', 'Position', [5 3 18 12])
scatter(results.plotting.modalities(:,8), results.plotting.modalities(:,2), ...
    'filled', '^'); % mean task data
ylim([-10 10]);
line('XData', [0 length(results.observers)], 'YData', [0, 0], 'LineStyle', '-', ...
    'LineWidth', 0.5, 'Color', 'k'); %midpoint
% Adding SD shaded area
ax = gca;
xVal = [ax.XLim(1):ax.XLim(end)];
shadedVal = zeros(length(xVal),1); %making the same length so can plot shaded error bar around 0
hold on
createShadedRegion(xVal, shadedVal, (shadedVal - SDpt5), (shadedVal + SDpt5),':','color', [0.1 0.1 0.1]);
hold on
createShadedRegion(xVal, shadedVal, (shadedVal - SD2), (shadedVal + SD2),':','color', [0.5 0.5 0.5]);
% Adding error bars to the mean data
hold on
errorbar(results.plotting.modalities(:,8), results.plotting.modalities(:,2), CIall, 'LineStyle', 'none',...
    'LineWidth', 0.7, 'Color', [0 0 0], 'CapSize', 0);
% Making it prettier
set(ax, 'FontSize', 10);
set(ax, 'XTick', results.observers);
xlabel('Observers'); ylabel('Bias (mm)');
pdfFileName = strcat('meanBias', '.pdf');
saveas(gcf, pdfFileName);

%% Making plot to show mean bias in all modalities
pdfFileName = strcat('biasModalities', '.pdf');
pngFileName = strcat('biasModalities', '.png');

figure('units', 'centimeters', 'Position', [5 3 18 12])
hold on
m1 = scatter(results.plotting.modalities(:,8), results.plotting.modalities(:,3), ...
    'filled', 'o', 'MarkerFaceColor', [0 0.7 0.2]); % landmark task data
set(m1, 'SizeData', 50);
hold on
m2 = scatter(results.plotting.modalities(:,8), results.plotting.modalities(:,4), ...
    'filled', 'o', 'MarkerFaceColor', [0.1 0 0.8]); % mlb task data
set(m2, 'SizeData', 50);
hold on
m3 = scatter(results.plotting.modalities(:,8), results.plotting.modalities(:,5), ...
    'filled', 'o', 'MarkerFaceColor', [0.4 0 0.4]); % trb task data
set(m3, 'SizeData', 50);
hold on
errorbar(results.plotting.modalities(:,8), results.plotting.modalities(:,2), CIall, 'LineStyle', 'none',...
    'LineWidth', 1, 'Color', [0 0 0], 'CapSize', 0);
hold on
m4 = scatter(results.plotting.modalities(:,8), results.plotting.modalities(:,2), ...
    'filled', '^', 'MarkerFaceColor', [0.6 0.6 0.9], 'MarkerEdgeColor', [0 0 0.1]); % mean task data
set(m4, 'SizeData', 60);
ylim([-15 15]);
line('XData', [0 length(results.observers)], 'YData', [0, 0], 'LineStyle', '-', ...
    'LineWidth', 0.5, 'Color', 'k'); %midpoint
% Adding SD shaded area
ax = gca;
xVal = [ax.XLim(1):ax.XLim(end)];
shadedVal = zeros(length(xVal),1); %making the same length so can plot shaded error bar around 0
hold on
createShadedRegion(xVal, shadedVal, (shadedVal - SDpt5), (shadedVal + SDpt5),':','color', [0.1 0.1 0]);
hold on
createShadedRegion(xVal, shadedVal, (shadedVal - SD2), (shadedVal + SD2),':','color', [0.7 0.7 0.5]);
% Adding error bars to the mean data
% Making it prettier
set(ax, 'FontSize', 10);
xLabels = num2str(results.plotting.modalities(:,1));
xticks(ax, 1:length(results.plotting.modalities(:,1)));
xticklabels(ax, xLabels);
%set(ax, 'XTick', results.observers);
xlabel('Observers'); ylabel('Bias (mm)');
mlgd = legend([m1 m2 m3 m4], 'Landmarks', 'MLB', 'TRB', 'Mean', [115 280 0.2 0.1]);
legend boxoff
mText = [mlgd, mlgd.ItemText]; set(mText, 'FontSize', 10);
% Adding text to define bias grouping
leftDim = [0.17 0.13 0.275 0.045];
rightDim = [0.765 0.13 0.12 0.045]; midDim = [0.45 0.13 0.31 0.045];
% annotation('textbox', leftDim, 'String', 'Left', 'FontSize', 8); %'Color', [0.2 0.2 0.2]);
% annotation('textbox', rightDim, 'String', 'Right', 'FontSize', 8); %'Color', [0.7 0.7 0.7]);
% annotation('textbox', midDim, 'String', 'Unspecified', 'FontSize', 8);  %'Color', [0.45 0.45 0.45]);
% saving as both pdf and png for ease
saveas(gcf, pdfFileName);
saveas(gcf, pngFileName);

%set(gca, 'Xtick', results.plotting.modalities(:,1));
% need to 'prettyfy', add correct labels, add group info, change axes, standard deviation, midline bar

%% Sessions plot
% Data matrix for sessions - all
% Landmarks
results.plotting.sessions.lm(:,1) = results.observers;
results.plotting.sessions.lm(:,2) = allData.modalities.data(:,1); %mean all data for lm
results.plotting.sessions.lm(:,3) = allData.sessions.lmPSE(:,1); %session 1
results.plotting.sessions.lm(:,4) = allData.sessions.lmPSE(:,2); %session 2
results.plotting.sessions.lm(:,5) = allData.sessions.lmPSE(:,3); %session 3
results.plotting.sessions.lm(:,6) = allData.sessions.lmPSE(:,4); %session 4
results.plotting.sessions.lm(:,7) = std(allData.sessions.lmPSE,0,2); %std of all sessions for lm
% Extracting the CIs from the SDs
lmCI = (results.plotting.sessions.lm(:,7)/sqrt(720))*1.97;
results.plotting.sessions.lm(:,8) = lmCI;
% Sorting the matrix by mean bias
results.plotting.sessions.lm = sortrows(results.plotting.sessions.lm, 2);
results.plotting.sessions.lm(:,9) = results.observers; %observers not sorted by mean bias for use with plotting - organisation of data
% standard deviation values for shading
lmSDpt5 = allData.means.lmPSE(2)*0.5;
lmSD2 = allData.means.lmPSE(2)*2;
lmSDall = results.plotting.sessions.lm(:,7);
lmCI = results.plotting.sessions.lm(:,8); %need to do this again because data has been sorted

% MLB
results.plotting.sessions.mlb(:,1) = results.observers;
results.plotting.sessions.mlb(:,2) = allData.modalities.data(:,2); %mean all data for mlb
results.plotting.sessions.mlb(:,3) = allData.sessions.mlb(:,1); %session 1
results.plotting.sessions.mlb(:,4) = allData.sessions.mlb(:,2); %session 2
results.plotting.sessions.mlb(:,5) = allData.sessions.mlb(:,3); %session 3
results.plotting.sessions.mlb(:,6) = allData.sessions.mlb(:,4); %session 4
results.plotting.sessions.mlb(:,7) = std(allData.sessions.mlb,0,2); %std of all sessions for mlb
% Extracting the CIs from the SDs
mlbCI = (results.plotting.sessions.mlb(:,7)/sqrt(360))*1.97;
results.plotting.sessions.mlb(:,8) = mlbCI;
% Sorting the matrix by mean bias
results.plotting.sessions.mlb = sortrows(results.plotting.sessions.mlb, 2);
results.plotting.sessions.mlb(:,9) = results.observers; %observers not sorted by mean bias for use with plotting - organisation of data
% standard deviation values for shading
mlbSDpt5 = allData.means.mlb(2)*0.5;
mlbSD2 = allData.means.mlb(2)*2;
mlbSDall = results.plotting.sessions.mlb(:,7);
mlbCI = results.plotting.sessions.mlb(:,8); %need to do this again because data has been sorted

% TRB
results.plotting.sessions.trb(:,1) = results.observers;
results.plotting.sessions.trb(:,2) = allData.modalities.data(:,3); %mean all data for trn
results.plotting.sessions.trb(:,3) = allData.sessions.trb(:,1); %session 1
results.plotting.sessions.trb(:,4) = allData.sessions.trb(:,2); %session 2
results.plotting.sessions.trb(:,5) = allData.sessions.trb(:,3); %session 3
results.plotting.sessions.trb(:,6) = allData.sessions.trb(:,4); %session 4
results.plotting.sessions.trb(:,7) = std(allData.sessions.trb,0,2); %std of all sessions for trb
% Extracting the CIs from the SDs
trbCI = (results.plotting.sessions.trb(:,7)/sqrt(216))*1.98;
results.plotting.sessions.trb(:,8) = trbCI;
% Sorting the matrix by mean bias
results.plotting.sessions.trb = sortrows(results.plotting.sessions.trb, 2);
results.plotting.sessions.trb(:,9) = results.observers; %observers not sorted by mean bias for use with plotting - organisation of data
% standard deviation values for shading
trbSDpt5 = allData.means.trb(2)*0.5;
trbSD2 = allData.means.trb(2)*2;
trbSDall = results.plotting.sessions.trb(:,7);
trbCI = results.plotting.sessions.trb(:,8); %need to do this again because data has been sorted

% All modalities - across sessions
results.plotting.sessions.all(:,1) = results.observers;
results.plotting.sessions.all(:,2) = allData.means.allPP(:,1); %mean all data
results.plotting.sessions.all(:,3) = allData.sessions.Session01all(:,1); %session 1
results.plotting.sessions.all(:,4) = allData.sessions.Session02all(:,1); %session 2
results.plotting.sessions.all(:,5) = allData.sessions.Session03all(:,1); %session 3
results.plotting.sessions.all(:,6) = allData.sessions.Session04all(:,1); %session 4
allData.sessions.allSessions = [allData.sessions.Session01all(:,1), allData.sessions.Session02all(:,1), ...
    allData.sessions.Session03all(:,1), allData.sessions.Session04all(:,1)];
results.plotting.sessions.all(:,7) = std(allData.sessions.allSessions,0,2); %std of all sessions for trb
% Extracting the CIs from the SDs
sessallCI = (results.plotting.sessions.all(:,7)/sqrt(19))*1.98;
results.plotting.sessions.all(:,8) = sessallCI;
% Sorting the matrix by mean bias
results.plotting.sessions.all = sortrows(results.plotting.sessions.all, 2);
results.plotting.sessions.all(:,9) = results.observers; %observers not sorted by mean bias for use with plotting - organisation of data
% standard deviation values for shading
allSDpt5 = std(allData.means.allPP(:,1))*0.5;
allSD2 = std(allData.means.allPP(:,1))*2;
allSDall = results.plotting.sessions.all(:,7);
allCI = results.plotting.sessions.all(:,8); %need to do this again because data has been sorted

%% Plotting landmarks all PP
pdfFileName = strcat('biasSessions_lm', '.pdf');
pngFileName = strcat('biasSessions_lm', '.png');

figure('units', 'centimeters', 'Position', [5 3 18 12])
lm1 = scatter(results.plotting.sessions.lm(:,9), results.plotting.sessions.lm(:,3), ...
    'filled', 'o', 'MarkerFaceColor', [0 0.2 0]); % landmark task data
set(lm1, 'SizeData', 50);
hold on
lm2 = scatter(results.plotting.sessions.lm(:,9), results.plotting.sessions.lm(:,4), ...
    'filled', 'o', 'MarkerFaceColor', [0 0.5 0]); % mlb task data
set(lm2, 'SizeData', 50);
hold on
lm3 = scatter(results.plotting.sessions.lm(:,9), results.plotting.sessions.lm(:,5), ...
    'filled', 'o', 'MarkerFaceColor', [0.1 0.7 0.1]); % trb task data
set(lm3, 'SizeData', 50);
hold on
lm4 = scatter(results.plotting.sessions.lm(:,9), results.plotting.sessions.lm(:,6), ...
    'filled', 'o', 'MarkerFaceColor', [0.4 1 0.4]); % trb task data
set(lm4, 'SizeData', 50);
% Adding error bars to the mean data
hold on
lm5 = scatter(results.plotting.sessions.lm(:,9), results.plotting.sessions.lm(:,2), ...
    'filled', '^', 'MarkerFaceColor', [0.8 1 0.6], 'MarkerEdgeColor', [0 0.2 0]); % mean task data
set(lm5, 'SizeData', 60);
hold on
errorbar(results.plotting.sessions.lm(:,9), results.plotting.sessions.lm(:,2), lmCI, 'LineStyle', 'none',...
    'LineWidth', 1, 'Color', [0 0 0], 'CapSize', 0);
ylim([-10 10]);
line('XData', [0 length(results.observers)], 'YData', [0, 0], 'LineStyle', '-', ...
    'LineWidth', 0.5, 'Color', 'k'); %midpoint
% Adding SD shaded area
ax = gca;
xVal = [ax.XLim(1):ax.XLim(end)];
shadedVal = zeros(length(xVal),1); %making the same length so can plot shaded error bar around 0
hold on
createShadedRegion(xVal, shadedVal, (shadedVal - lmSDpt5), (shadedVal + lmSDpt5),':','color', [0.1 0.1 0]);
hold on
createShadedRegion(xVal, shadedVal, (shadedVal - lmSD2), (shadedVal + lmSD2),':','color', [0.7 0.7 0.5]);
% Making it prettier
set(ax, 'FontSize', 14);
xLabels = num2str(results.plotting.sessions.lm(:,1));
xticks(ax, 1:length(results.plotting.sessions.lm(:,1)));
xticklabels(ax, xLabels);
xlabel('Observers'); ylabel('Bias (mm)');
lmlgd = legend([lm1, lm2, lm3, lm4, lm5], '1', '2', '3', '4', 'Mean', [105 262 0.2 0.1]);
lmText = [lmlgd, lmlgd.ItemText]; set(lmText, 'FontSize', 12);
legend boxoff
% Adding text to define bias grouping
leftDim = [0.17 0.13 0.235 0.045];
rightDim = [0.575 0.13 0.29 0.045]; midDim = [0.41 0.13 0.16 0.045];
% annotation('textbox', leftDim, 'String', 'Left', 'FontSize', 8); %'Color', [0.2 0.2 0.2]);
% annotation('textbox', rightDim, 'String', 'Right', 'FontSize', 8); %'Color', [0.7 0.7 0.7]);
% annotation('textbox', midDim, 'String', 'Unspecified', 'FontSize', 8);  %'Color', [0.45 0.45 0.45]);
% % saving as both pdf and png for ease
saveas(gcf, pdfFileName);
saveas(gcf, pngFileName);

%% Plotting mlb all PP
pdfFileName = strcat('biasSessions_mlb', '.pdf');
pngFileName = strcat('biasSessions_mlb', '.png');

figure('units', 'centimeters', 'Position', [5 3 18 12])
mlb1 = scatter(results.plotting.sessions.mlb(:,9), results.plotting.sessions.mlb(:,3), ...
    'filled', 'o', 'MarkerFaceColor', [0 0 0.2]); % landmark task data
set(mlb1, 'SizeData', 50);
hold on
mlb2 = scatter(results.plotting.sessions.mlb(:,9), results.plotting.sessions.mlb(:,4), ...
    'filled', 'o', 'MarkerFaceColor', [0 0 0.5]); % mlb task data
set(mlb2, 'SizeData', 50);
hold on
mlb3 = scatter(results.plotting.sessions.mlb(:,9), results.plotting.sessions.mlb(:,5), ...
    'filled', 'o', 'MarkerFaceColor', [0.1 0.1 0.7]); % trb task data
set(mlb3, 'SizeData', 50);
hold on
mlb4 = scatter(results.plotting.sessions.lm(:,9), results.plotting.sessions.mlb(:,6), ...
    'filled', 'o', 'MarkerFaceColor', [0.5 0.5 1]); % 
set(mlb4, 'SizeData', 50);
hold on
mlb5 = scatter(results.plotting.sessions.mlb(:,9), results.plotting.sessions.mlb(:,2), ...
    'filled', '^', 'MarkerFaceColor', [0.8 0.8 1], 'MarkerEdgeColor', [0 0 0.1]); % mean task data
set(mlb5, 'SizeData', 60);
ylim([-25 25]);
line('XData', [0 length(results.observers)], 'YData', [0, 0], 'LineStyle', '-', ...
    'LineWidth', 0.5, 'Color', 'k'); %midpoint
% Adding error bars to the mean data
hold on
errorbar(results.plotting.sessions.mlb(:,9), results.plotting.sessions.mlb(:,2), mlbCI, 'LineStyle', 'none',...
    'LineWidth', 1, 'Color', [0 0 0], 'CapSize', 0);
% Adding SD shaded area
ax = gca;
xVal = [ax.XLim(1):ax.XLim(end)];
shadedVal = zeros(length(xVal),1); %making the same length so can plot shaded error bar around 0
hold on
createShadedRegion(xVal, shadedVal, (shadedVal - mlbSDpt5), (shadedVal + mlbSDpt5),':','color', [0.1 0.1 0]);
hold on
createShadedRegion(xVal, shadedVal, (shadedVal - mlbSD2), (shadedVal + mlbSD2),':','color', [0.7 0.7 0.5]);
% Making it prettier
set(ax, 'FontSize', 14);
xLabels = num2str(results.plotting.sessions.mlb(:,1));
xticks(ax, 1:length(results.plotting.sessions.mlb(:,1)));
xticklabels(ax, xLabels);
xlabel('Observers'); ylabel('Bias (mm)');
mlblgd = legend([mlb1, mlb2, mlb3, mlb4, mlb5], '1', '2', '3', '4', 'Mean', [105 262 0.2 0.1]);
mlbText = [mlblgd, mlblgd.ItemText]; set(mlbText, 'FontSize', 12);
legend boxoff
% Adding text to define bias grouping
leftDim = [0.17 0.13 0.39 0.045];
rightDim = [0.82 0.13 0.06 0.045]; midDim = [0.565 0.13 0.25 0.045];
% annotation('textbox', leftDim, 'String', 'Left', 'FontSize', 8); %'Color', [0.2 0.2 0.2]);
% annotation('textbox', rightDim, 'String', 'Right', 'FontSize', 8); %'Color', [0.7 0.7 0.7]);
% annotation('textbox', midDim, 'String', 'Unspecified', 'FontSize', 8);  %'Color', [0.45 0.45 0.45]);
% % saving as both pdf and png for ease
saveas(gcf, pdfFileName);
saveas(gcf, pngFileName);

%% Plotting trb all PP
pdfFileName = strcat('biasSessions_trb', '.pdf');
pngFileName = strcat('biasSessions_trb', '.png');

figure('units', 'centimeters', 'Position', [5 3 18 12])
trb1 = scatter(results.plotting.sessions.trb(:,9), results.plotting.sessions.trb(:,3), ...
    'filled', 'o', 'MarkerFaceColor', [0.3 0.2 0.1]); % landmark task data
set(trb1, 'SizeData', 50);
hold on
trb2 = scatter(results.plotting.sessions.trb(:,9), results.plotting.sessions.trb(:,4), ...
    'filled', 'o', 'MarkerFaceColor', [0.6 0.3 0.2]); % mlb task data
set(trb2, 'SizeData', 50);
hold on
trb3 = scatter(results.plotting.sessions.trb(:,9), results.plotting.sessions.trb(:,5), ...
    'filled', 'o', 'MarkerFaceColor', [0.85 0.5 0.3]); % trb task data
set(trb3, 'SizeData', 50);
hold on
trb4 = scatter(results.plotting.sessions.trb(:,9), results.plotting.sessions.trb(:,6), ...
    'filled', 'o', 'MarkerFaceColor', [1 0.7 0.4]); % trb task data
set(trb4, 'SizeData', 50);
hold on
trb5 = scatter(results.plotting.sessions.trb(:,9), results.plotting.sessions.trb(:,2), ...
    'filled', '^', 'MarkerFaceColor', [1 0.9 0.7], 'MarkerEdgeColor', [0.2 0.1 0]); % mean task data
set(trb5, 'SizeData', 60);
ylim([-25 25]);
% Adding error bars to the mean data
hold on
errorbar(results.plotting.sessions.trb(:,9), results.plotting.sessions.trb(:,2), trbCI, 'LineStyle', 'none',...
    'LineWidth', 1, 'Color', [0 0 0], 'CapSize', 0);
line('XData', [0 length(results.observers)], 'YData', [0, 0], 'LineStyle', '-', ...
    'LineWidth', 0.5, 'Color', 'k'); %midpoint
% Adding SD shaded area
ax = gca;
xVal = [ax.XLim(1):ax.XLim(end)];
shadedVal = zeros(length(xVal),1); %making the same length so can plot shaded error bar around 0
hold on
createShadedRegion(xVal, shadedVal, (shadedVal - trbSDpt5), (shadedVal + trbSDpt5),':','color', [0.1 0.1 0]);
hold on
createShadedRegion(xVal, shadedVal, (shadedVal - trbSD2), (shadedVal + trbSD2),':','color', [0.7 0.7 0.5]);
% Making it prettier
set(ax, 'FontSize', 14);
xLabels = num2str(results.plotting.sessions.trb(:,1));
xticks(ax, 1:length(results.plotting.sessions.trb(:,1)));
xticklabels(ax, xLabels);
xlabel('Observers'); ylabel('Bias (mm)');
trblgd = legend([trb1, trb2, trb3, trb4, trb5], '1', '2', '3', '4', 'Mean', [105 262 0.2 0.1]);
trbText = [trblgd, trblgd.ItemText]; set(trbText, 'FontSize', 12);
legend boxoff
% Adding text to define bias grouping
leftDim = [0.17 0.13 0.27 0.045];
rightDim = [0.68 0.13 0.195 0.045]; midDim = [0.445 0.13 0.23 0.045];
% annotation('textbox', leftDim, 'String', 'Left', 'FontSize', 8); %'Color', [0.2 0.2 0.2]);
% annotation('textbox', rightDim, 'String', 'Right', 'FontSize', 8); %'Color', [0.7 0.7 0.7]);
% annotation('textbox', midDim, 'String', 'Unspecified', 'FontSize', 8);  %'Color', [0.45 0.45 0.45]);
% saving as both pdf and png for ease
saveas(gcf, pdfFileName);
saveas(gcf, pngFileName);

%% Plotting mean sessions all PP
pdfFileName = strcat('biasSessions_all', '.pdf');
pngFileName = strcat('biasSessions_all', '.png');

figure('units', 'centimeters', 'Position', [5 3 18 12])
all1 = scatter(results.plotting.sessions.all(:,9), results.plotting.sessions.all(:,3), ...
    'filled', 'o', 'MarkerFaceColor', [0.2 0 0.2]); % landmark task data
set(all1, 'SizeData', 50);
hold on
all2 = scatter(results.plotting.sessions.all(:,9), results.plotting.sessions.all(:,4), ...
    'filled', 'o', 'MarkerFaceColor', [0.5 0 0.5]); % mlb task data
set(all2, 'SizeData', 50);
hold on
all3 = scatter(results.plotting.sessions.all(:,9), results.plotting.sessions.all(:,5), ...
    'filled', 'o', 'MarkerFaceColor', [0.7 0.2 0.8]); % trb task data
set(all3, 'SizeData', 50);
hold on
all4 = scatter(results.plotting.sessions.all(:,9), results.plotting.sessions.all(:,6), ...
    'filled', 'o', 'MarkerFaceColor', [1 0.4 1]); % trb task data
set(all4, 'SizeData', 50);
hold on
all5 = scatter(results.plotting.sessions.all(:,9), results.plotting.sessions.all(:,2), ...
    'filled', '^', 'MarkerFaceColor', [1 0.8 1], 'MarkerEdgeColor', [0.1 0 0.1]); % mean task data
set(all5, 'SizeData', 60);
ylim([-15 15]);
% Adding error bars to the mean data
hold on
errorbar(results.plotting.sessions.all(:,9), results.plotting.sessions.all(:,2), allCI, 'LineStyle', 'none',...
    'LineWidth', 1, 'Color', [0 0 0], 'CapSize', 0);
line('XData', [0 length(results.observers)], 'YData', [0, 0], 'LineStyle', '-', ...
    'LineWidth', 0.5, 'Color', 'k'); %midpoint
% Adding SD shaded area
ax = gca;
xVal = [ax.XLim(1):ax.XLim(end)];
shadedVal = zeros(length(xVal),1); %making the same length so can plot shaded error bar around 0
hold on
createShadedRegion(xVal, shadedVal, (shadedVal - allSDpt5), (shadedVal + allSDpt5),':','color', [0.1 0.1 0]);
hold on
createShadedRegion(xVal, shadedVal, (shadedVal - allSD2), (shadedVal + allSD2),':','color', [0.7 0.7 0.5]);
% Making it prettier
set(ax, 'FontSize', 10);
xLabels = num2str(results.plotting.sessions.all(:,1));
xticks(ax, 1:length(results.plotting.sessions.all(:,1)));
xticklabels(ax, xLabels);
xlabel('Observers'); ylabel('Bias (mm)');
allLgd = legend([all1, all2, all3, all4, all5], '1', '2', '3', '4', 'Mean', [105 275 0.2 0.1]);
allText = [allLgd, allLgd.ItemText]; set(allText, 'FontSize', 10);
legend boxoff
% Adding text to define bias grouping
leftDim = [0.17 0.13 0.27 0.045];
rightDim = [0.68 0.13 0.195 0.045]; midDim = [0.445 0.13 0.23 0.045];
% annotation('textbox', leftDim, 'String', 'Left', 'FontSize', 8); %'Color', [0.2 0.2 0.2]);
% annotation('textbox', rightDim, 'String', 'Right', 'FontSize', 8); %'Color', [0.7 0.7 0.7]);
% annotation('textbox', midDim, 'String', 'Unspecified', 'FontSize', 8);  %'Color', [0.45 0.45 0.45]);
% saving as both pdf and png for ease
saveas(gcf, pdfFileName);
saveas(gcf, pngFileName);

%% Plotting modality mean across all sessions
% Going to try this 2 ways
% First - bias y axis, session xaxis, modalities separate
xValues = 1:(length(nSessions));
sessSD1 = std(allData.sessions.allSessions(:,1));
sessSD2 = std(allData.sessions.allSessions(:,2));
sessSD3 = std(allData.sessions.allSessions(:,3));
sessSD4 = std(allData.sessions.allSessions(:,4));
sessCI1 = (sessSD1/sqrt(length(results.observers)))*2.11;
sessCI2 = (sessSD2/sqrt(length(results.observers)))*2.11;
sessCI3 = (sessSD3/sqrt(length(results.observers)))*2.11;
sessCI4 = (sessSD4/sqrt(length(results.observers)))*2.11;
errorBars = [sessCI1, sessCI2, sessCI3, sessCI4];
allData.sessions.means.all = mean([allData.sessions.means.lmPSE; allData.sessions.means.mlb; ...
    allData.sessions.means.trb]);

pdfFileName = strcat('biasSessions_bySess', '.pdf');
pngFileName = strcat('biasSessions_bySess', '.png');

figure('units', 'centimeters', 'Position', [5 3 12 10])
dataMean = line('XData', [0 length(xValues)], 'YData', [allData.means.tot(1), allData.means.tot(1)],...
    'LineStyle', '-', 'LineWidth', 2, 'Color', [0.7 0.7 0.7]); % mean for all data
hold on
errorbar(xValues, allData.sessions.means.all, errorBars, 'LineStyle', 'none', 'LineWidth', 1.5,...
    'Color', [0.4 0.4 0.4], 'CapSize', 0);
hold on
s1 = scatter(xValues, allData.sessions.means.lmPSE, 'o', 'MarkerFaceColor', [0 0.8 0.1], ...
    'MarkerEdgeColor', [0 0.8 0.1]); %lm data
set(s1, 'SizeData', 60);
hold on
s2 = scatter(xValues, allData.sessions.means.mlb, 'filled', 'square', 'MarkerFaceColor', [0.2 0.1 0.7], ...
    'MarkerEdgeColor', [0.2 0.1 0.7]);
set(s2, 'SizeData', 60);
hold on
s3 = scatter(xValues, allData.sessions.means.trb, 'filled', '^', 'MarkerFaceColor', [0.5 0 0.5], ...
    'MarkerEdgeColor', [0.5 0 0.5]);
set(s3, 'SizeData', 60);
midpoint = line('XData', [0 length(xValues)], 'YData', [0, 0], 'LineStyle', '--', ...
    'LineWidth', 0.5, 'Color', [0 0 0]); %midpoint
ylim([-6 6]);
ax = gca;
set(ax, 'FontSize', 10, 'xtick', [1 2 3 4]);
xLabels = {'1', '2', '3', '4'};
xticks(ax, [1 2 3 4]);
xticklabels(ax, xLabels);
ylabel('Bias (mm)');
xlabel('Sessions');
f3lgd = legend([s1 s2 s3], 'Landmarks', 'MLB', 'TRB', [88 233 0.1 0.3]);
f3Text = [f3lgd, f3lgd.ItemText]; set(f3Text, 'FontSize', 10);
legend boxoff
saveas(gcf, pdfFileName);
saveas(gcf, pngFileName);


% Second - bias y axis, modality xaxis, session separate
% Values
xValues = 1:3;
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));
    allData.modalities.means.(sprintf('%s', session))... 
        = [allData.sessions.means.lmPSE(i), allData.sessions.means.mlb(i), allData.sessions.means.trb(i)];
end
allData.modalities.means.all = mean([allData.modalities.means.Session01; allData.modalities.means.Session02;...
    allData.modalities.means.Session03; allData.modalities.means.Session04]);
allData.modalities.means.sdall = std([allData.modalities.means.Session01; allData.modalities.means.Session02;...
    allData.modalities.means.Session03; allData.modalities.means.Session04]);
allData.modalities.means.CIall = (allData.modalities.means.sdall/sqrt(length(results.observers)))*2.11;
    
pdfFileName = strcat('biasSessions_byMod', '.pdf');
pngFileName = strcat('biasSessions_byMod', '.png');

figureS = figure('units', 'centimeters', 'Position', [5 3 12 10]);
dataMean = line('XData', [0 length(xValues)], 'YData', [allData.means.tot(1), allData.means.tot(1)], 'LineStyle', '-', ...
    'LineWidth', 2, 'Color', [0.7 0.7 0.7;]); % mean for all data
hold on
s1 = scatter(xValues, allData.modalities.means.Session01, 'diamond', 'MarkerFaceColor', [0.1 0 0.1], ...
    'MarkerEdgeColor', [0.1 0 0.1]); %lm data
set(s1, 'SizeData', 60);
hold on
errorbar(xValues, allData.modalities.means.all, allData.modalities.means.CIall, 'LineStyle', 'none', 'LineWidth', 1.5,...
    'Color', [0.4 0.4 0.4], 'CapSize', 0);
hold on
s2 = scatter(xValues, allData.modalities.means.Session02, 'filled', 'o', 'MarkerFaceColor', [0.4 0 0.4], ...
    'MarkerEdgeColor', [0.3 0 0.3]);
set(s2, 'SizeData', 60);
hold on
s3 = scatter(xValues, allData.modalities.means.Session03, '^', 'MarkerFaceColor', [0.7 0 0.7], ...
    'MarkerEdgeColor', [0.5 0 0.5]);
set(s3, 'SizeData', 60);
hold on
s4 = scatter(xValues, allData.modalities.means.Session04, 'square', 'MarkerFaceColor', [0.9 0.2 0.9], ...
    'MarkerEdgeColor', [0.7 0 0.7]);
set(s4, 'SizeData', 60);
midpoint = line('XData', [0 length(xValues)], 'YData', [0, 0], 'LineStyle', '--', ...
    'LineWidth', 0.5, 'Color', [0 0 0]); %midpoint
ylim([-6 6]);
ax = gca;
set(ax, 'FontSize', 10);
xLabels = {'Landmarks', 'MLB', 'TRB'};
xticks(ax, [1 2 3]);
xticklabels(ax, xLabels);
%xlabel('Modalities'); 
ylabel('Bias (mm)'); xlabel('Task');
f4lgd = legend([s1, s2, s3, s4], '1', '2', '3', '4', [70 225 0.1 0.3]);
legend boxoff
f4Text = [f4lgd, f4lgd.ItemText]; set(f4Text, 'FontSize', 10);
saveas(gcf, pdfFileName);
saveas(gcf, pngFileName);

%% Cronbach's alpha and other stats...
type = 'C-k'; %type of ICC used - 2-k, 2-way fixed effects
% Analysing reliability across sessions
[results.sessions.lm.r, results.sessions.lm.bound(1), results.sessions.lm.bound(2), results.sessions.lm.F,...
    results.sessions.lm.df(1), results.sessions.lm.df(2), results.sessions.lm.p] = ICC(allData.sessions.lmPSE, type);
[results.sessions.mlb.r, results.sessions.mlb.bound(1), results.sessions.mlb.bound(2), results.sessions.mlb.F,...
    results.sessions.mlb.df(1), results.sessions.mlb.df(2), results.sessions.mlb.p] = ICC(allData.sessions.mlb, type);
[results.sessions.trb.r, results.sessions.trb.bound(1), results.sessions.trb.bound(2), results.sessions.trb.F,...
    results.sessions.trb.df(1), results.sessions.trb.df(2), results.sessions.trb.p] = ICC(allData.sessions.trb, type);
% Reliability across all sessions
[results.sessions.all.r, results.sessions.all.bound(1), results.sessions.all.bound(2), results.sessions.all.F,...
    results.sessions.all.df(1), results.sessions.all.df(2), results.sessions.all.p] = ICC(allData.sessions.allSessions, type);

% Reliability across modalities
[results.modalities.all.r, results.modalities.all.bound(1), results.modalities.all. bound(2), results.modalities.all.F,...
    results.modalities.all.df(1), results.modalities.all.df(2), results.modalities.all.p] = ICC(allData.modalities.data, type);

% One-sample t-tests
[results.modalities.lmT.h results.modalities.lmT.p results.modalities.lmT.ci ...
    results.modalities.lmT.stats] = ttest(allData.modalities.data(:,1)); %landmark one-sample
[results.modalities.mlbT.h results.modalities.mlbT.p results.modalities.mlbT.ci ...
    results.modalities.mlbT.stats] = ttest(allData.modalities.data(:,2)); %mlb one-sample
[results.modalities.trbT.h results.modalities.trbT.p results.modalities.trbT.ci ...
    results.modalities.trbT.stats] = ttest(allData.modalities.data(:,3)); %trb one-sample
[results.all.h results.all.p results.all.ci results.all.stats] = ...
    ttest(allData.means.allPP(:,1));

% bonferroni corrected
results.modalities.lmT.Bp = results.modalities.lmT.p*3;
results.modalities.mlbT.Bp = results.modalities.mlbT.p*3;
results.modalities.trbT.Bp = results.modalities.trbT.p*3;

%% Creating covariance matrices
% Landmarks
% Getting  matrix for all session data
lmcovar = cov(allData.sessions.lmPSE); %creating covariance matrix
lmcovarCorr = corrcov(lmcovar);
lmcovarScaled = mat2gray(lmcovarCorr, [1 0]);
imwrite(lmcovarScaled, 'lmcovar', 'PNG'); %writing to an image
lmcovarImg = imread('lmcovar'); %reading the image back in, so it's in img format

lmTextStr = num2str(lmcovarCorr(:), '%0.2f'); %creating a text string for placing over the image
lmTextStr = strtrim(cellstr(lmTextStr)); %remove any funny spacings
[x, y] = meshgrid(1:4); %coordinates for the numbers

% Plotting matrix, plot all together
pngFileName = 'CovarMatrix.png';
figure('pos',[150 150 1200 300])
subplot(1,3,1)
imagesc(lmcovarScaled);
colormap gray
C = colormap; L = size(C,1);
Gs = round(interp1(linspace(min(lmcovarScaled(:)), max(lmcovarScaled(:)),L),1:L,lmcovarScaled));
H = reshape(C(Gs,:), [size(Gs) 3]);
image(H)
hold on
lmStrings = text(x(:), y(:), lmTextStr(:), 'HorizontalAlignment', 'center');
ax = gca;
set(ax, 'FontSize', 12);
xLabels = {'Sess 1', 'Sess 2', 'Sess 3', 'Sess 4'};
yLabels = {'Sess 1', 'Sess 2', 'Sess 3', 'Sess 4'};
xticks(ax, [1 2 3 4]); yticks(ax, [1 2 3 4]);
xticklabels(ax, xLabels); yticklabels(ax, yLabels);
ytickangle(90)
title('Landmarks')

%% MLB
mlbcovar = cov(allData.sessions.mlb); %creating covariance matrix
mlbcovarCorr = corrcov(mlbcovar);
mlbcovarScaled = mat2gray(mlbcovarCorr, [1 0]);
imwrite(mlbcovarScaled, 'mlbcovar', 'PNG'); %writing to an image
mlbcovarImg = imread('mlbcovar'); %reading the image back in, so it's in img format

mlbTextStr = num2str(mlbcovarCorr(:), '%0.2f'); %creating a text string for placing over the image
mlbTextStr = strtrim(cellstr(mlbTextStr)); %remove any funny spacings

subplot(1,3,2)
imagesc(mlbcovarScaled);
colormap gray
C = colormap; L = size(C,1);
Gs = round(interp1(linspace(min(mlbcovarScaled(:)), max(mlbcovarScaled(:)),L),1:L,mlbcovarScaled));
H = reshape(C(Gs,:), [size(Gs) 3]);
image(H)
hold on
mlbStrings = text(x(:), y(:), mlbTextStr(:), 'HorizontalAlignment', 'center');
ax = gca;
set(ax, 'FontSize', 12);
xLabels = {'Sess 1', 'Sess 2', 'Sess 3', 'Sess 4'};
yLabels = {' ', '', '', ''};
xticks(ax, [1 2 3 4]); yticks(ax, [1 2 3 4]);
xticklabels(ax, xLabels); yticklabels(ax, yLabels);
ytickangle(90)
title('Manual line bisection')

%% TRB
trbcovar = cov(allData.sessions.trb); %creating covariance matrix
trbcovarCorr = corrcov(trbcovar);
trbcovarScaled = mat2gray(trbcovarCorr, [1 0]);
imwrite(trbcovarScaled, 'trbcovar', 'PNG'); %writing to an image
trbcovarImg = imread('trbcovar'); %reading the image back in, so it's in img format

trbTextStr = num2str(trbcovarCorr(:), '%0.2f'); %creating a text string for placing over the image
trbTextStr = strtrim(cellstr(trbTextStr)); %remove any funny spacings

subplot(1,3,3)
imagesc(trbcovarScaled);
colormap gray
C = colormap; L = size(C,1);
Gs = round(interp1(linspace(min(trbcovarScaled(:)), max(trbcovarScaled(:)),L),1:L,trbcovarScaled));
H = reshape(C(Gs,:), [size(Gs) 3]);
image(H)
hold on
trbStrings = text(x(:), y(:), trbTextStr(:), 'HorizontalAlignment', 'center');
ax = gca;
set(ax, 'FontSize', 12);
xLabels = {'Sess 1', 'Sess 2', 'Sess 3', 'Sess 4'};
yLabels = {' ', '', '', ''};
xticks(ax, [1 2 3 4]); yticks(ax, [1 2 3 4]);
xticklabels(ax, xLabels); yticklabels(ax, yLabels);
ytickangle(90)
title('Tactile rod bisection')
saveas(gcf, pngFileName);

figure(4)
ax = gca;
set(ax, 'FontSize', 14);
colormap gray
colorbar('eastoutside', 'Direction', 'reverse')
saveas(figure(4), 'covarColourbar.png');


%% Modalities
modcovar = cov(allData.modalities.data);
modcovarCorr = corrcov(modcovar);

modTextStr = num2str(modcovarCorr(:), '%0.2f'); %creating a text string for placing over the image
modTextStr = strtrim(cellstr(modTextStr)); %remove any funny spacings
[xM, yM] = meshgrid(1:3);

pngFileName = 'CovarMatrix_mod.png';
figure()
imagesc(modcovarCorr)
colormap(flipud(gray))
hold on
modStrings = text(xM(:), yM(:), modTextStr(:), 'HorizontalAlignment', 'center',...
    'FontSize', 14);
ax = gca;
set(ax, 'FontSize', 14);
xLabels = {'LM', 'MLB', 'TRB'};
yLabels = {'LM', 'MLB', 'TRB'};
xticks(ax, [1 2 3]); yticks(ax, [1 2 3]);
xticklabels(ax, xLabels); yticklabels(ax, yLabels);
ytickangle(90)
saveas(gcf, pngFileName);

%% save and close
close all
cd(dirAnaAll)
save(matfilename, 'allData', 'results');
