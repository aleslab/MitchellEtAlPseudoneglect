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
    modlm = mean(allData.sessions.lmPSE(p,:));
    modlmSD = std(allData.sessions.lmPSE(p,:));
    allData.modalities.data(p,1) = modlm; allData.modalities.sds(p,1) = modlmSD;
    % Manual line bisection
    modmlb = mean(allData.sessions.mlb(p,:));
    modmlbSD = std(allData.sessions.mlbProp(p,:));
    allData.modalities.data(p,2) = modmlb; allData.modalities.sds(p,2) = modmlbSD;
    % Tactile rod bisection
    modtrb = mean(allData.sessions.trb(p,:));
    modtrbSD = std(allData.sessions.trb(p,:));
    allData.modalities.data(p,3) = modtrb; allData.modalities.sds(p,3) = modtrbSD;
    
    
    %% Mean for all PP
    % column 1 means, column 2 SDs
    allData.means.allPP(p,1) = mean([modlm, modmlb, modtrb]);
    allData.means.allPP(p,2) = std([modlm, modmlb, modtrb]);
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
    
    % Modalities across sessions
    allData.modalities.data(remove(r),:) = NaN;
    allData.modalities.sds(remove(r),:) = NaN;
    
    % All means
    allData.means.allPP(remove(r),:) = NaN;
end
% Removing outlier from matrix - modalities and sessions
allData.sessions.lmPSE = allData.sessions.lmPSE(find(allData.sessions.lmPSE(:,1) > 0.00001 | ...
    allData.sessions.lmPSE(:,1) < 0.00001),:);
allData.sessions.CIhigh= allData.sessions.CIhigh(find(allData.sessions.CIhigh(:,1) > 0.00001 | ...
    allData.sessions.CIhigh(:,1) < 0.00001),:);
allData.sessions.CIlow= allData.sessions.CIlow(find(allData.sessions.CIlow(:,1) > 0.00001 | ...
    allData.sessions.CIlow(:,1) < 0.00001),:);
allData.sessions.lmStd= allData.sessions.lmStd(find(allData.sessions.lmStd(:,1) > 0.00001 | ...
    allData.sessions.lmStd(:,1) < 0.00001),:);
allData.sessions.mlb= allData.sessions.mlb(find(allData.sessions.mlb(:,1) > 0.00001 | ...
    allData.sessions.mlb(:,1) < 0.00001),:);
allData.sessions.mlbStd= allData.sessions.mlbStd(find(allData.sessions.mlbStd(:,1) > 0.00001 | ...
    allData.sessions.mlbStd(:,1) < 0.00001),:);
allData.sessions.trb= allData.sessions.trb(find(allData.sessions.trb(:,1) > 0.00001 | ...
    allData.sessions.trb(:,1) < 0.00001),:);
allData.sessions.trbStd= allData.sessions.trbStd(find(allData.sessions.trbStd(:,1) > 0.00001 | ...
    allData.sessions.trbStd(:,1) < 0.00001),:);

% Modalities across sessions
allData.modalities.data = allData.modalities.data(find(allData.modalities.data(:,1) > 0.00001 | ...
    allData.modalities.data(:,1) < 0.00001),:);
allData.modalities.sds = allData.modalities.sds(find(allData.modalities.sds(:,1) > 0.00001 | ...
    allData.modalities.sds(:,1) < 0.00001),:);

% All means
allData.means.allPP = allData.means.allPP(find(allData.means.allPP(:,1) > 0.00001 | ...
    allData.means.allPP(:,1) < 0.00001),:);

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
% Sorting the matrix by mean bias
results.plotting.modalities = sortrows(results.plotting.modalities, 2);
results.plotting.modalities(:,7) = results.observers; %observers not sorted by mean bias for use with plotting - organisation of data

% standard deviation values for shading
SDpt5 = allData.means.tot(2)*0.5;
SD2 = allData.means.tot(2)*2;
SDall = results.plotting.modalities(:,6);

cd(dirAnaAll);
% Making initial mean error (all modalities) plot
figure('units', 'centimeters', 'Position', [5 3 18 12])
scatter(results.plotting.modalities(:,7), results.plotting.modalities(:,2), ...
    'filled', '^'); % mean task data
ylim([-3 3]);
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
errorbar(results.plotting.modalities(:,7), results.plotting.modalities(:,2), SDall, 'LineStyle', 'none',...
    'LineWidth', 0.7, 'Color', [0 0 0], 'CapSize', 0);
% Making it prettier
set(ax, 'FontSize', 12);
set(ax, 'XTick', results.observers);
xlabel('Observers'); ylabel('Bias (mm)');
pdfFileName = strcat('meanBias', '.pdf');
saveas(gcf, pdfFileName);

%% Making plot to show bias in all modalities
pdfFileName = strcat('biasModalities', '.pdf');
pngFileName = strcat('biasModalities', '.png');

figure('units', 'centimeters', 'Position', [5 3 18 12])
hold on
scatter(results.plotting.modalities(:,7), results.plotting.modalities(:,3), ...
    'filled', 'o', 'MarkerFaceColor', [0 0 0.4]); % landmark task data
hold on
scatter(results.plotting.modalities(:,7), results.plotting.modalities(:,4), ...
    'filled', 'o', 'MarkerFaceColor', [0 0 0.8]); % mlb task data
hold on
scatter(results.plotting.modalities(:,7), results.plotting.modalities(:,5), ...
    'filled', 'o', 'MarkerFaceColor', [0.6 0.6 1]); % trb task data
hold on
scatter(results.plotting.modalities(:,7), results.plotting.modalities(:,2), ...
    'filled', '^', 'MarkerEdgeColor', [0.3 0 0.3]); % mean task data
ylim([-3 3]);
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
errorbar(results.plotting.modalities(:,7), results.plotting.modalities(:,2), SDall, 'LineStyle', 'none',...
    'LineWidth', 0.7, 'Color', [0 0 0], 'CapSize', 0);
% Making it prettier
set(ax, 'FontSize', 12);
set(ax, 'XTick', results.observers);
xlabel('Observers'); ylabel('Bias (mm)');
legend('Landmarks', 'MLB', 'TRB', 'Mean', [120 280 0.2 0.1]);
% Adding text to define bias grouping
leftDim = [0.17 0.13 0.25 0.045];
rightDim = [0.77 0.13 0.13 0.045]; midDim = [0.43 0.13 0.33 0.045];
annotation('textbox', leftDim, 'String', 'Left', 'FontSize', 8); %'Color', [0.2 0.2 0.2]);
annotation('textbox', rightDim, 'String', 'Right', 'FontSize', 8); %'Color', [0.7 0.7 0.7]);
annotation('textbox', midDim, 'String', 'Unspecified', 'FontSize', 8);  %'Color', [0.45 0.45 0.45]);
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
% Sorting the matrix by mean bias
results.plotting.sessions.lm = sortrows(results.plotting.sessions.lm, 2);
results.plotting.sessions.lm(:,8) = results.observers; %observers not sorted by mean bias for use with plotting - organisation of data
% standard deviation values for shading
lmSDpt5 = allData.means.lmPSE(2)*0.5;
lmSD2 = allData.means.lmPSE(2)*2;
lmSDall = results.plotting.sessions.lm(:,7);

% MLB
results.plotting.sessions.mlb(:,1) = results.observers;
results.plotting.sessions.mlb(:,2) = allData.modalities.data(:,2); %mean all data for mlb
results.plotting.sessions.mlb(:,3) = allData.sessions.mlb(:,1); %session 1
results.plotting.sessions.mlb(:,4) = allData.sessions.mlb(:,2); %session 2
results.plotting.sessions.mlb(:,5) = allData.sessions.mlb(:,3); %session 3
results.plotting.sessions.mlb(:,6) = allData.sessions.mlb(:,4); %session 4
results.plotting.sessions.mlb(:,7) = std(allData.sessions.mlb,0,2); %std of all sessions for lm
% Sorting the matrix by mean bias
results.plotting.sessions.mlb = sortrows(results.plotting.sessions.mlb, 2);
results.plotting.sessions.mlb(:,8) = results.observers; %observers not sorted by mean bias for use with plotting - organisation of data
% standard deviation values for shading
mlbSDpt5 = allData.means.mlb(2)*0.5;
mlbSD2 = allData.means.mlb(2)*2;
mlbSDall = results.plotting.sessions.mlb(:,7);

% MLB
results.plotting.sessions.trb(:,1) = results.observers;
results.plotting.sessions.trb(:,2) = allData.modalities.data(:,3); %mean all data for trn
results.plotting.sessions.trb(:,3) = allData.sessions.trb(:,1); %session 1
results.plotting.sessions.trb(:,4) = allData.sessions.trb(:,2); %session 2
results.plotting.sessions.trb(:,5) = allData.sessions.trb(:,3); %session 3
results.plotting.sessions.trb(:,6) = allData.sessions.trb(:,4); %session 4
results.plotting.sessions.trb(:,7) = std(allData.sessions.trb,0,2); %std of all sessions for lm
% Sorting the matrix by mean bias
results.plotting.sessions.trb = sortrows(results.plotting.sessions.trb, 2);
results.plotting.sessions.trb(:,8) = results.observers; %observers not sorted by mean bias for use with plotting - organisation of data
% standard deviation values for shading
trbSDpt5 = allData.means.trb(2)*0.5;
trbSD2 = allData.means.trb(2)*2;
trbSDall = results.plotting.sessions.trb(:,7);


%% Plotting landmarks all PP
pdfFileName = strcat('biasSessions_lm', '.pdf');
pngFileName = strcat('biasSessions_lm', '.png');

figure('units', 'centimeters', 'Position', [5 3 18 12])
hold on
scatter(results.plotting.sessions.lm(:,8), results.plotting.sessions.lm(:,3), ...
    'filled', 'o', 'MarkerFaceColor', [0 0.2 0]); % landmark task data
hold on
scatter(results.plotting.sessions.lm(:,8), results.plotting.sessions.lm(:,4), ...
    'filled', 'o', 'MarkerFaceColor', [0 0.5 0]); % mlb task data
hold on
scatter(results.plotting.sessions.lm(:,8), results.plotting.sessions.lm(:,5), ...
    'filled', 'o', 'MarkerFaceColor', [0.1 0.7 0.1]); % trb task data
hold on
scatter(results.plotting.sessions.lm(:,8), results.plotting.sessions.lm(:,6), ...
    'filled', 'o', 'MarkerFaceColor', [0.5 1 0.5]); % trb task data
hold on
scatter(results.plotting.sessions.lm(:,8), results.plotting.sessions.lm(:,2), ...
    'filled', '^', 'MarkerEdgeColor', [0.1 0 0]); % mean task data
ylim([-4 4]);
line('XData', [0 length(results.observers)], 'YData', [0, 0], 'LineStyle', '-', ...
    'LineWidth', 0.5, 'Color', 'k'); %midpoint
% Adding SD shaded area
ax = gca;
xVal = [ax.XLim(1):ax.XLim(end)];
shadedVal = zeros(length(xVal),1); %making the same length so can plot shaded error bar around 0
hold on
createShadedRegion(xVal, shadedVal, (shadedVal - lmSDpt5), (shadedVal + lmSDpt5),':','color', [0.1 0.1 0.1]);
hold on
createShadedRegion(xVal, shadedVal, (shadedVal - lmSD2), (shadedVal + lmSD2),':','color', [0.5 0.5 0.5]);
% Adding error bars to the mean data
hold on
errorbar(results.plotting.sessions.lm(:,8), results.plotting.sessions.lm(:,2), lmSDall, 'LineStyle', 'none',...
    'LineWidth', 0.7, 'Color', [0 0 0], 'CapSize', 0);
% Making it prettier
set(ax, 'FontSize', 12);
set(ax, 'XTick', results.observers);
xlabel('Observers'); ylabel('Bias (mm)');
legend('1', '2', '3', '4', 'Mean', [110 270 0.2 0.1]);
% Adding text to define bias grouping
leftDim = [0.17 0.13 0.225 0.045];
rightDim = [0.64 0.13 0.26 0.045]; midDim = [0.40 0.13 0.235 0.045];
annotation('textbox', leftDim, 'String', 'Left', 'FontSize', 8); %'Color', [0.2 0.2 0.2]);
annotation('textbox', rightDim, 'String', 'Right', 'FontSize', 8); %'Color', [0.7 0.7 0.7]);
annotation('textbox', midDim, 'String', 'Unspecified', 'FontSize', 8);  %'Color', [0.45 0.45 0.45]);
% saving as both pdf and png for ease
saveas(gcf, pdfFileName);
saveas(gcf, pngFileName);

%% Plotting modality mean across all sessions
% Going to try this 2 ways
% First - bias y axis, session xaxis, modalities separate
xValues = 1:length(nSessions);
figure()
scatter(xValues, allData.sessions.means.lmPSE, 'o', 'MarkerFaceColor', [0 0.8 0.1]); %lm data
hold on
scatter(xValues, allData.sessions.means.mlb, 'filled', 'o', 'MarkerFaceColor', [0.2 0.1 0.7]);
hold on
scatter(xValues, allData.sessions.means.trb, 'filled', '^', 'MarkerFaceColor', [0.5 0 0.5]);
midpoint = line('XData', [0 length(xValues)], 'YData', [0, 0], 'LineStyle', '-', ...
    'LineWidth', 0.5, 'Color', 'k'); %midpoint
hold on
dataMean = line('XData', [0 length(xValues)], 'YData', [allData.means.tot(1), allData.means.tot(1)], 'LineStyle', '-', ...
    'LineWidth', 0.8, 'Color', [0.5 0.5 0.5]); % mean for all data
ylim([-1 1]);

% Second - bias y axis, modality xaxis, session separate

%% Cronbach's alpha and other stats...
type = 'C-k'; %type of ICC used - 2-k, 2-way fixed effects
% Analysing reliability across sessions
[results.sessions.lm.r, results.sessions.lm.bound(1), results.sessions.lm.bound(2), results.sessions.lm.F,...
    results.sessions.lm.df(1), results.sessions.lm.df(2), results.sessions.lm.p] = ICC(allData.sessions.lmPSE, type);
[results.sessions.mlb.r, results.sessions.mlb.bound(1), results.sessions.mlb.bound(2), results.sessions.mlb.F,...
    results.sessions.mlb.df(1), results.sessions.mlb.df(2), results.sessions.mlb.p] = ICC(allData.sessions.mlb, type);
[results.sessions.trb.r, results.sessions.trb.bound(1), results.sessions.trb.bound(2), results.sessions.trb.F,...
    results.sessions.trb.df(1), results.sessions.trb.df(2), results.sessions.trb.p] = ICC(allData.sessions.trb, type);

% Reliability across modalities
[results.modalities.all.r, results.modalities.all.bound(1), results.modalities.all. bound(2), results.modalities.all.F,...
    results.modalities.all.df(1), results.modalities.all.df(2), results.modalities.all.p] = ICC(allData.modalities.data, type);

%% save and close
close all
cd(dirAnaAll)
save(matfilename, 'allData', 'results');
