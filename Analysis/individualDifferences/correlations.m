%% AG.Mitchell - 19.06.19
%% Exploratory analysis
% Investigating correlations between standard bisection and 2AFC tasks to
% see whether data from one can predict the other in the same participant
% Landmark/line bisection correlated with landmark 2AFC, tactile rod
% bisection correlated with tactile rod 2AFC

%% Directories and variables
clear all
close all

dataDir = 'M:\Alex_Files\Experiments\Bias\Analysis';
reliabilityFile = 'ReliabilityAnalysis.mat';
responseFile = 'responseAnalysis.mat';
cd(dataDir)

% loading reliability data
load(reliabilityFile)
allData_rel = allData;
results_rel = results;
% loading response data
load(responseFile)
allData_res = allData;
results_res = results;

% Removing outliers from 2AFC that are present in the standard bisection
% Sorting back into pp order
lm_resData1 = sortrows(results_res.plotting.shifts.lm2,1);
tr_resData1 = sortrows(results_res.plotting.shifts.tr2,1);
% outlier removal
remove = find(outlierSum > 0.1); %identifies paticipants that need removing
remove = sort(remove); %sorting so in participant order
for r = 1:length(remove)
    lm_resData1(remove(r),:) = NaN;
    tr_resData1(remove(r),:) = NaN;
end

% Getting data sets without NaN values
lm_resData = lm_resData1(isfinite(lm_resData1(:,1)),:);
tr_resData = tr_resData1(isfinite(tr_resData1(:,1)),:);

% Getting data sets for bisection tasks (and sorting back into PP order)
lm_relData = sortrows(results_rel.plotting.sessions.lm,1);
mlb_relData = sortrows(results_rel.plotting.sessions.mlb,1);
trb_relData = sortrows(results_rel.plotting.sessions.trb,1);

%% Landmark vs. LM2AFC
% Plotting first
% Need to calculate polyfits for line of best fit
regressionCoeff_P = polyfit(lm_relData(:,2), lm_resData(:,2),1);
lmX_plot = [min(lm_relData(:,2)); max(lm_relData(:,2))]; %making data plotable
lmY_plot_P = polyval(regressionCoeff_P, lmX_plot);

lmMark = 50;
% Take the mean values across all testing sessions
figure()
lmS = scatter(lm_relData(:,2), lm_resData(:,2), lmMark, 'filled', 'MarkerEdgeColor', [0.3 0.1 0.4],...
    'MarkerFaceColor', [0.7 0.4 0.8]);
grid on
% axes
ax = gca;
ylim([0.3 0.6]); xlim([-5 4])
%set(ax, 'FontSize', 11, 'xtick', -5:1:5, 'ytick', 0:0.9);
% y axes labels = 2dp
ytix = get(ax,'ytick')';
xtix = get(ax,'xtick')';
set(ax,'yticklabel',num2str(ytix,'%.2f'))
% drawing lines at 0 point
line('XData', [0 0], 'YData', [ytix(1) ytix(end)], 'LineStyle', '-', ...
    'LineWidth', 1, 'Color', [0.4 0.4 0.4]); %x-axis
line('XData', [xtix(1) xtix(end)], 'YData', [0.5 0.5], 'LineStyle', '-', ...
    'LineWidth', 1, 'Color', [0.4 0.4 0.4]); %x-axis
% drawing line of best fit
hold on
lmP = plot(lmX_plot, lmY_plot_P, 'Color', [0.3, 0.1, 0.4], 'LineWidth', 2); %polyfit
% labelling axes
lmT = title('Landmark bisection vs. 2AFC'); set(lmT, 'FontSize', 12);
lmAxT(1) = xlabel('Landmark bisection error (mm)');
lmAxT(2) = ylabel('Landmark 2AFC proportion RHS perceived longer');

% Correlation second

