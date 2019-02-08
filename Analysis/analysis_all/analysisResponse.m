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
        allData.(sprintf('%s', session)).lm2.allPP(:,(p)) = (lm2psych.(sprintf('%s', session))(:,2))/100; %calculating proportion
        % tactile data
        allData.(sprintf('%s', session)).tr2.asymmetry = tr2psych.(sprintf('%s', session))(:,1);
        allData.(sprintf('%s', session)).tr2.allPP(:,(p)) = (tr2psych.(sprintf('%s', session))(:,2))/100;
        
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
    % tactile
    allData.allSessions.tr2.asymmetry = tr2psych.allSessions(:,1);
    allData.allSessions.tr2.allPP(:,(p)) = (tr2psych.allSessions(:,2))/100;
    
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
    results.lm2.binomial{p,:} = true(allData.allSessions.lm2.shiftGroup(1,p)<results.lm2.binomialStat(1)); %landmarks
    results.tr2.binomial{p,:} = true(allData.allSessions.tr2.shiftGroup(1,p)<results.tr2.binomialStat(1));
end

%% Plotting mean data
% Make a plot that runs from -10 to 0 asymmetry as both lines shifted by
% the same amount

%% Landmarks
cd(dirAnaAll) %directory where all analysis stored

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
%% Plot Binomial