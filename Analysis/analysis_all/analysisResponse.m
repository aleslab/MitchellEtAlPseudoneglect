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
dirBias = 'M:\Alex_Files\Experiments\Bias\'; %subject to change depending on where you analyse
nParticipants = [1:19,21:24,26:30];
%nParticipants = 3; %for testing
nSessions = [1:4];
allData = struct;

% Getting outlier information
dirAnaAll = [dirBias filesep 'Analysis']; %directory for all analysis - here is where data should be saved from this file
cd(dirAnaAll)
dataFilename = ('ReliabilityAnalysis.mat');
load(dataFilename)

for p = 1:length(nParticipants)
    ppID = sprintf('P%0*d',2,nParticipants(p)); %for use when navigating files
    % Variables 
    nSessions = 1:4;
    visualFilename = sprintf('%s_visualanalysis2AFC.mat', ppID);
    tactileFilename = sprintf('%s_tactileanalysis.mat', ppID);
    matFilename = ('responseAnalysis.mat');
    % Directory
    dirPP = [dirBias 'Data' filesep ppID]; %participant directory
    dirAna = [dirPP filesep 'Analysis' filesep];
    dirVis = [dirAna 'Visual' filesep];
    dirTact = [dirAna 'Tactile' filesep];
    
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
        allData.(sprintf('%s', session)).lm2.allPP(:,(p)) = (lm2psych.(sprintf('%s', session))(:,2))/100; %calculating proportion
        allData.(sprintf('%s', session)).lm2.SD(p) = std(allData.(sprintf('%s', session)).lm2.allPP(1:6,p)); %sd of all 6 shifts (-10 to 0)
        % tactile data
        allData.(sprintf('%s', session)).tr2.asymmetry = tr2psych.(sprintf('%s', session))(:,1);
        allData.(sprintf('%s', session)).tr2.allPP(:,(p)) = (tr2psych.(sprintf('%s', session))(:,2))/100;
        allData.(sprintf('%s', session)).tr2.SD(p) = std(allData.(sprintf('%s', session)).tr2.allPP(1:6,p));
        
        % average data into the 3 largest and 3 smallest stimulus asymmetries
        % landmarks
        largeShift = allData.(sprintf('%s', session)).lm2.allPP(1:3,(p));
        allData.(sprintf('%s', session)).lm2.shiftGroup(1,p) = mean(largeShift); %proportions for -10 to -6mm shift
        smallShift = allData.(sprintf('%s', session)).lm2.allPP(4:6,(p));
        allData.(sprintf('%s', session)).lm2.shiftGroup(2,p) = mean(smallShift); %proportions for -4 to 0mm shift  
        % tactile
        largeShift = allData.(sprintf('%s', session)).tr2.allPP(1:3,(p));
        allData.(sprintf('%s', session)).tr2.shiftGroup(1,p) = mean(largeShift); %proportions for -10 to -6mm shift
        smallShift = allData.(sprintf('%s', session)).tr2.allPP(4:6,(p));
        allData.(sprintf('%s', session)).tr2.shiftGroup(2,p) = mean(smallShift); %proportions for -4 to 0mm shift  
    end
    
    % Matrices outside sessions- average across all sessions for all PP
    allData.(sprintf('%s', session)).names{p} = ppID; %adding participant ID to name
    % landmarks
    allData.allSessions.lm2.asymmetry = lm2psych.allSessions(:,1); 
    allData.allSessions.lm2.allPP(:,(p)) = (lm2psych.allSessions(:,2))/100; %calculating proportions
    allData.allSessions.lm2.SD(p) = std(allData.allSessions.lm2.allPP(1:6,p));
    % tactile
    allData.allSessions.tr2.asymmetry = tr2psych.allSessions(:,1);
    allData.allSessions.tr2.allPP(:,(p)) = (tr2psych.allSessions(:,2))/100;
    allData.allSessions.tr2.SD(p) = std(allData.allSessions.tr2.allPP(1:6,p));
    
    % average data into the 3 largest and 3 smallest stimulus asymmetries
    % landmarks
    largeShift = allData.allSessions.lm2.allPP(1:3,(p));
    allData.allSessions.lm2.shiftGroup(1,p) = mean(largeShift); %proportions for -10 to -6mm shift
    smallShift = allData.allSessions.lm2.allPP(4:6,(p));
    allData.allSessions.lm2.shiftGroup(2,p) = mean(smallShift); %proportions for -4 to 0mm shift  
    % tactile
    largeShift = allData.allSessions.tr2.allPP(1:3,(p));
    allData.allSessions.tr2.shiftGroup(1,p) = mean(largeShift); %proportions for -10 to -6mm shift
    smallShift = allData.allSessions.tr2.allPP(4:6,(p));
    allData.allSessions.tr2.shiftGroup(2,p) = mean(smallShift); %proportions for -4 to 0mm shift  
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
    
    % big and small shifts
    % landmarks
    allData.(sprintf('%s', session)).lm2.shiftGroupAv(:,1) = ...
        mean(allData.(sprintf('%s', session)).lm2.shiftGroup, 2);
    allData.(sprintf('%s', session)).lm2.shiftGroupAv(:,2) = ...
        std(allData.(sprintf('%s', session)).lm2.shiftGroup, 0,2);
    % tactile
    allData.(sprintf('%s', session)).tr2.shiftGroupAv(:,1) = ...
        mean(allData.(sprintf('%s', session)).tr2.shiftGroup, 2);
    allData.(sprintf('%s', session)).tr2.shiftGroupAv(:,2) = ...
        std(allData.(sprintf('%s', session)).tr2.shiftGroup, 0,2);
end
% across all sessions
allData.allSessions.lm2.av(:,1) = mean(allData.allSessions.lm2.allPP, 2); %mean 
allData.allSessions.lm2.av(:,2) = std(allData.allSessions.lm2.allPP, 0,2); %std

allData.allSessions.tr2.av(:,1) = mean(allData.allSessions.tr2.allPP, 2); %mean 
allData.allSessions.tr2.av(:,2) = std(allData.allSessions.tr2.allPP, 0,2); %std

% big and small shifts
allData.allSessions.lm2.shiftGroupAv(:,1) = mean(allData.allSessions.lm2.shiftGroup, 2);
allData.allSessions.lm2.shiftGroupAv(:,2) = std(allData.allSessions.lm2.shiftGroup, 0,2);

allData.allSessions.tr2.shiftGroupAv(:,1) = mean(allData.allSessions.tr2.shiftGroup, 2);
allData.allSessions.tr2.shiftGroupAv(:,2) = std(allData.allSessions.tr2.shiftGroup, 0,2);

%% Binomial tests 
% For average across largest three shifts 
% For all participants

LMflips = 360; %number of trials for each participant for largest 3 shifts, to use in binomial
TRflips = 108;
alpha = 0.5; %testing against 50%

results.lm2.binomialStat(1) = binoinv(0.025, LMflips, alpha)/LMflips; %low
results.lm2.binomialStat(2) = binoinv(0.975, LMflips, alpha)/LMflips; %high
results.tr2.binomialStat(1) = binoinv(0.025, TRflips, alpha)/TRflips; %low
results.tr2.binomialStat(2) = binoinv(0.975, TRflips, alpha)/TRflips; %high

% Calculate whether each participant is below the boundary set by binomial
% data
for p = 1:length(nParticipants)
    ppID = sprintf('P%0*d',2,nParticipants(p)); %for use when navigating files
    
    % determining whether each participant is below the lower binomial
    % bounadry (therefore significantly different from 50%)
    % Test for left bias (large shift smaller than low binomial)
    results.lm2.binomial.logL{p,:} = true(allData.allSessions.lm2.shiftGroup(1,p)<results.lm2.binomialStat(1)); %landmarks
    results.tr2.binomial.logL{p,:} = true(allData.allSessions.tr2.shiftGroup(1,p)<results.tr2.binomialStat(1));
    % Test for right bias (large shift larger than high binomial)
    results.lm2.binomial.logR{p,:} = true(allData.allSessions.lm2.shiftGroup(1,p)>results.lm2.binomialStat(2)); %landmarks
    results.tr2.binomial.logR{p,:} = true(allData.allSessions.tr2.shiftGroup(1,p)>results.tr2.binomialStat(2));
end

% getting proportion of participants
% landmarks - left shift
results.lm2.binomial.trueL = double(cell2mat(results.lm2.binomial.logL));
results.lm2.binomial.proportionSigLeft = sum(results.lm2.binomial.trueL)/length(nParticipants);
% right shift
results.lm2.binomial.trueR = double(cell2mat(results.lm2.binomial.logR));
results.lm2.binomial.proportionSigRight = sum(results.lm2.binomial.trueR)/length(nParticipants);

% tactile - left shift
results.tr2.binomial.trueL = double(cell2mat(results.tr2.binomial.logL));
results.tr2.binomial.proportionSigLeft = sum(results.tr2.binomial.trueL)/length(nParticipants);
% right shift
results.tr2.binomial.trueR = double(cell2mat(results.tr2.binomial.logR));
results.tr2.binomial.proportionSigRight = sum(results.tr2.binomial.trueR)/length(nParticipants); 

%% Plotting mean data
% Make a plot that runs from -10 to 0 asymmetry as both lines shifted by
% the same amount

%% Landmarks
cd(dirAnaAll) %directory where all analysis stored
% for shaded regions
lm2low = 0.5-results.lm2.binomialStat(1); lm2high = results.lm2.binomialStat(2)-0.5;
tr2low = 0.5-results.tr2.binomialStat(1); tr2high = results.tr2.binomialStat(2)-0.5;

pdfFileName = strcat('lm2_responseAveragePP_sessions', '.pdf');
asym = allData.allSessions.lm2.asymmetry(1:6);
figure()
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));
    % plotting
    M = plot(asym, allData.(sprintf('%s', session)).lm2.av(1:6,1), 'LineWidth', 1.5);
    hold on
end
ax = gca;
set(ax, 'Xtick', asym); ylim([0 0.7]);
set(ax, 'FontSize', 11)
hold on %drawing a line at 50%
midpoint = line('XData', [-10 0], 'YData', [0.5, 0.5], 'LineStyle', '--', ...
    'LineWidth', 0.5, 'Color', [0 0 0]); %midpoint
xlabel('Stimulus asymmetry (mm)'); ylabel('Proportion right-shifted line perceived as longer');
% Include shaded area with binomial data
shadedVal = zeros(length(asym),1)+0.5; %making the same length so can plot shaded error bar around 0
hold on
createShadedRegion(asym, shadedVal, (shadedVal - lm2low),...
    (shadedVal + lm2high),':','color', [0.5 0.5 0.5]);
saveas(gcf, pdfFileName);

% Average across all sessions
% Getting the CI of the data
lmCI = (allData.allSessions.lm2.av(:,2)/sqrt(length(nParticipants)))*2.11;
trCI = (allData.allSessions.tr2.av(:,2)/sqrt(length(nParticipants)))*2.11;

% Plotting figure
figure()
pdfFileName = strcat('lm2_responseAveragePP', '.pdf');
plot(asym, allData.allSessions.lm2.av(1:6,1), 'LineWidth', 1);
hold on %plotting error bars
E = errorbar(asym, allData.allSessions.lm2.av(1:6,1), lmCI(1:6), 'Color', ...
    [0.3 0.3 0.3], 'LineWidth', 1.5);
ax = gca;
set(ax, 'Xtick', asym); ylim([0 0.7]);
set(ax, 'FontSize', 11);
hold on %drawing a line at 50%
midpoint = line('XData', [-10 0], 'YData', [0.5, 0.5], 'LineStyle', '--', ...
    'LineWidth', 0.5, 'Color', [0 0 0]); %midpoint
shadedVal = zeros(length(asym),1)+0.5; %making the same length so can plot shaded error bar around 0
hold on
createShadedRegion(asym, shadedVal, (shadedVal - lm2low),...
    (shadedVal + lm2high),':','color', [0.5 0.5 0.5]);
xlabel('Stimulus asymmetry (mm)'); ylabel('Proportion right-shifted line perceived as longer');
saveas(gcf, pdfFileName);

%% Tactile rod
pdfFileName = strcat('tr2_responseAveragePP_sessions', '.pdf');
asym = allData.allSessions.tr2.asymmetry(1:6);
figure()
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));
    % plotting
    M = plot(asym, allData.(sprintf('%s', session)).tr2.av(1:6,1), 'LineWidth', 1.5);
    hold on
end
ax = gca;
set(ax, 'Xtick', asym); ylim([0 0.7]);
set(ax, 'FontSize', 11);
hold on %drawing a line at 50%
midpoint = line('XData', [-10 0], 'YData', [0.5, 0.5], 'LineStyle', '--', ...
    'LineWidth', 0.5, 'Color', [0 0 0]); %midpoint
shadedVal = zeros(length(asym),1)+0.5; %making the same length so can plot shaded error bar around 0
hold on
createShadedRegion(asym, shadedVal, (shadedVal - tr2low),...
    (shadedVal + tr2high),':','color', [0.5 0.5 0.5]);
xlabel('Stimulus asymmetry (mm)'); ylabel('Proportion right-shifted line perceived as longer');
saveas(gcf, pdfFileName);

% Average across all sessions
figure()
pdfFileName = strcat('tr2_responseAveragePP', '.pdf');
plot(asym, allData.allSessions.tr2.av(1:6,1), 'LineWidth', 1.5)
hold on
E = errorbar(asym, allData.allSessions.tr2.av(1:6,1), trCI(1:6), 'Color', ...
    [0.3 0.3 0.3], 'LineWidth', 1.5);
ax = gca;
set(ax, 'Xtick', asym); ylim([0 0.7]);
set(ax, 'FontSize', 11);
hold on %drawing a line at 50%
midpoint = line('XData', [-10 0], 'YData', [0.5, 0.5], 'LineStyle', '--', ...
    'LineWidth', 0.5, 'Color', [0 0 0]); %midpoint
shadedVal = zeros(length(asym),1)+0.5; %making the same length so can plot shaded error bar around 0
hold on
createShadedRegion(asym, shadedVal, (shadedVal - tr2low),...
    (shadedVal + tr2high),':','color', [0.5 0.5 0.5]);
xlabel('Stimulus asymmetry (mm)'); ylabel('Proportion right-shifted line perceived as longer');
saveas(gcf, pdfFileName);

%% Out of interest - plot the grouped shifts for both tasks
data = [allData.allSessions.lm2.shiftGroupAv(1,1), allData.allSessions.lm2.shiftGroupAv(2,1); ...
    allData.allSessions.tr2.shiftGroupAv(1,1), allData.allSessions.tr2.shiftGroupAv(2,1)];
error = [allData.allSessions.lm2.shiftGroupAv(1,2), allData.allSessions.lm2.shiftGroupAv(2,2); ...
    allData.allSessions.tr2.shiftGroupAv(1,2), allData.allSessions.tr2.shiftGroupAv(2,2)];

figure()
B = bar(data, 'BarWidth', 0.8);
B(1,1).FaceColor = [0.4 0.4 0.4]; B(1,2).FaceColor = [0.4 0.3 0.7];
% setting error bars
ax = gca;
xticks(ax, [1,2]); xticklabels(ax, {'Landmarks 2AFC', 'Tactile Rod 2AFC'});
ylabel('Proportion right-shifted line perceived longer');
nGroups = size(data,1); %for error bar information
nBars = size(data,2);
% group width for each bar in the group
groupWidth = min(0.8, nBars/(nBars + 1.5));
hold on %adding error bar
for i = 1:nBars
    x = (1:nGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*nBars);
    EB(:,i) = errorbar(x, data(:,i), error(:,i), 'LineWidth', 1);
end
EB(1,1).Color = [0.2 0.2 0.2]; EB(1,2).Color = [0.2 0.1 0.5];
legend('Large shift', 'Small shift');

%% Making Dakin plot
%% Landmarks
% Plotting the difference between average shifts
results.plotting.shifts.lm2(:,1) = nParticipants; %participant info
results.plotting.shifts.lm2(:,2) = mean(allData.allSessions.lm2.shiftGroup); %average proportion
results.plotting.shifts.lm2(:,3) = allData.allSessions.lm2.shiftGroup(1,:); %high shift
results.plotting.shifts.lm2(:,4) = allData.allSessions.lm2.shiftGroup(2,:); %low shift
results.plotting.shifts.lm2(:,5) = allData.allSessions.lm2.SD; %std of all shifts
% Extracting CIs from the SDs
lm2CI = results.plotting.shifts.lm2(:,5)/sqrt(length(nParticipants))*2.10;
results.plotting.shifts.lm2(:,6) = lm2CI;
% organising participants by mean proportion
results.plotting.shifts.lm2 = sortrows(results.plotting.shifts.lm2, 2);
results.plotting.shifts.lm2(:,7) = 1:22; %observers not sorted by bias for plotting
% standard deviation values for shading
lm2SDpt5 = std(results.plotting.shifts.lm2(:,2))*0.5;
lm2SD2 = std(results.plotting.shifts.lm2(:,2))*2;
lm2SDall = results.plotting.shifts.lm2(:,5);
lm2CI = results.plotting.shifts.lm2(:,6); %need to do this again because data has been sorted

% Making plot
pdfFileName = strcat('Landmarks2_groups', '.pdf');
figure('units', 'centimeters', 'Position', [5 3 18 12])
hold on
lm1 = scatter(results.plotting.shifts.lm2(:,7), results.plotting.shifts.lm2(:,3), ...
    'filled', 'o', 'MarkerFaceColor', [0.7 0.5 0.7], 'MarkerEdgeColor', [0.1 0 0.1]); % landmark task data
set(lm1, 'SizeData', 50);
hold on
lm2 = scatter(results.plotting.shifts.lm2(:,7), results.plotting.shifts.lm2(:,4), ...
    'filled', 'o', 'MarkerFaceColor', [0.3 0.1 0.3], 'MarkerEdgeColor', [0.3 0.2 0.3]); % mlb task data
set(lm2, 'SizeData', 50);
hold on %adding error bars
lmE = errorbar(results.plotting.shifts.lm2(:,7), results.plotting.shifts.lm2(:,2), lm2CI,  'LineStyle', 'none',...
    'LineWidth', 1, 'Color', [0 0 0], 'CapSize', 0);
hold on
tr3 = scatter(results.plotting.shifts.lm2(:,7), results.plotting.shifts.lm2(:,2), ...
    'filled', '^', 'MarkerFaceColor', [0.4 0 0.4], 'MarkerEdgeColor', [0.2 0 0.2]); % mean of task data
set(tr3, 'SizeData', 50);
hold on %drawing a line at 0.5
line('XData', [0 length(nParticipants)], 'YData', [0.5, 0.5], 'LineStyle', '-', ...
    'LineWidth', 0.5, 'Color', 'k'); %midpoint
% adding the fiddly bits
ax = gca;
xVal = [ax.XLim(1):ax.XLim(end)];
shadedVal = zeros(length(xVal),1)+0.5; %making the same length so can plot shaded error bar around 0
hold on %shaded region of binomial analysis
createShadedRegion(xVal, shadedVal, (shadedVal - lm2low), (shadedVal + lm2high),':','color', [0.7 0.7 0.7]);
% Adding SD shaded area
ylim([0.2 0.7]);
% Making it prettier
set(ax, 'FontSize', 11);
xLabels = num2str(results.plotting.shifts.lm2(:,1));
xticks(ax, 1:length(results.plotting.shifts.lm2(:,1)));
xticklabels(ax, xLabels); 
tix = get(gca, 'ytick')';
set(gca, 'yticklabel', num2str(tix, '%.2f')); %setting all to 1dp
% labels and legends
xlabel('Observers'); ylabel('Proportion right-shifted longer');
lmlgd = legend([lm1, lm2, tr3], 'Large shift', 'Small shift', 'Mean', [110 290 0.2 0.1]);
lmText = [lmlgd, lmlgd.ItemText]; set(lmText, 'FontSize', 10);
legend boxoff
saveas(gcf, pdfFileName);

%% Tactile
% Plotting the difference between average shifts
results.plotting.shifts.tr2(:,1) = nParticipants; %participant info
results.plotting.shifts.tr2(:,2) = mean(allData.allSessions.tr2.shiftGroup); %average proportion
results.plotting.shifts.tr2(:,3) = allData.allSessions.tr2.shiftGroup(1,:); %high shift
results.plotting.shifts.tr2(:,4) = allData.allSessions.tr2.shiftGroup(2,:); %low shift
results.plotting.shifts.tr2(:,5) = allData.allSessions.tr2.SD; %std of all shifts
% Extracting CIs from the SDs
tr2CI = results.plotting.shifts.tr2(:,5)/sqrt(length(nParticipants))*2.10;
results.plotting.shifts.tr2(:,6) = tr2CI;
% organising participants by mean proportion
results.plotting.shifts.tr2 = sortrows(results.plotting.shifts.tr2, 2);
results.plotting.shifts.tr2(:,7) = 1:22; %observers not sorted by bias for plotting
% standard deviation values for shading
tr2SDpt5 = std(results.plotting.shifts.tr2(:,2))*0.5;
tr2SD2 = std(results.plotting.shifts.tr2(:,2))*2;
tr2SDall = results.plotting.shifts.tr2(:,5);
tr2CI = results.plotting.shifts.tr2(:,6); %need to do this again because data has been sorted

% Making plot
pdfFileName = strcat('Tactile2_groups', '.pdf');
figure('units', 'centimeters', 'Position', [5 3 18 12])
hold on
tr1 = scatter(results.plotting.shifts.tr2(:,7), results.plotting.shifts.tr2(:,3), ...
    'filled', 'o', 'MarkerFaceColor', [0.2 0.9 0.2], 'MarkerEdgeColor', [0 0.5 0]); % landmark task data
set(tr1, 'SizeData', 50);
hold on
tr2 = scatter(results.plotting.shifts.tr2(:,7), results.plotting.shifts.tr2(:,4), ...
    'filled', 'o', 'MarkerFaceColor', [0 0.4 0], 'MarkerEdgeColor', [0 0.2 0]); % mlb task data
set(tr2, 'SizeData', 50);
hold on %adding error bars
trE = errorbar(results.plotting.shifts.tr2(:,7), results.plotting.shifts.tr2(:,2), tr2CI,  'LineStyle', 'none',...
    'LineWidth', 1, 'Color', [0 0 0], 'CapSize', 0);
hold on
tr3 = scatter(results.plotting.shifts.tr2(:,7), results.plotting.shifts.tr2(:,2), ...
    'filled', '^', 'MarkerFaceColor', [0 0.2 0], 'MarkerEdgeColor', [0 0.1 0]); % mean of task data
set(tr3, 'SizeData', 50);
hold on %drawing a line at 0.5
line('XData', [0 length(nParticipants)], 'YData', [0.5, 0.5], 'LineStyle', '-', ...
    'LineWidth', 0.5, 'Color', 'k'); %midpoint
% adding the fiddly bits
ax = gca;
xVal = [ax.XLim(1):ax.XLim(end)];
shadedVal = zeros(length(xVal),1)+0.5; %making the same length so can plot shaded error bar around 0
hold on %shaded region of binomial analysis
createShadedRegion(xVal, shadedVal, (shadedVal - tr2low), (shadedVal + tr2high),':','color', [0.7 0.7 0.7]);
% Adding SD shaded area
ylim([0.2 0.7]);
% Making it prettier
set(ax, 'FontSize', 11);
xLabels = num2str(results.plotting.shifts.tr2(:,1));
xticks(ax, 1:length(results.plotting.shifts.tr2(:,1)));
xticklabels(ax, xLabels); 
tix = get(gca, 'ytick')';
set(gca, 'yticklabel', num2str(tix, '%.2f')); %setting all to 1dp
% labels and legends
xlabel('Observers'); ylabel('Proportion right-shifted longer');
trlgd = legend([tr1, tr2, tr3], 'Large shift', 'Small shift', 'Mean', [110 290 0.2 0.1]);
trText = [trlgd, trlgd.ItemText]; set(trText, 'FontSize', 10);
legend boxoff
saveas(gcf, pdfFileName);

%% T-test
% Calculating difference between large and small shift group for each task
% Within participants
% Landmarks
[h, p, ci, stat] = ttest(results.plotting.shifts.lm2(:,3), results.plotting.shifts.lm2(:,4));
results.lm2.ttest.h = h; results.lm2.ttest.p = p; results.lm2.ttest.ci = ci;
results.lm2.ttest.stat = stat; %adding to results struct

% Tactile
[h, p, ci, stat] = ttest(results.plotting.shifts.tr2(:,3), results.plotting.shifts.tr2(:,4));
results.tr2.ttest.h = h; results.tr2.ttest.p = p; results.tr2.ttest.ci = ci;
results.tr2.ttest.stat = stat; %adding to results struct

save(matFilename, 'allData', 'results');
close all