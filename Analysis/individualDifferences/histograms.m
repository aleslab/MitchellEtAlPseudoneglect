%% AG. Mitchell - 14.06.19
%% Plotting individual participant data into a histogram
% Do different tasks show different spread of data?
% Look acros individual tasks and also individual testing sessions

% Take data from the reliabilityAnalysis.mat file (extracted using
% ReliabilityAnalysis.m script)

%% Data
% pathways
expDir = 'M:\Alex_Files\Experiments\Bias\Analysis\';
dataDir = [expDir 'analysis_all']; %original data set
resDir = [expDir 'individual_differences']; %for saving plots/individual diff analysis

% getting data
cd(dataDir)
load('reliabilityAnalysis.mat');

% variables
nSessions = 1:4;
nParticipants = 1:length(results.plotting.modalities(:,1));

% Data in easy-to-read variables
landmark.mean = results.plotting.sessions.lm(:,2);
lineBisection.mean = results.plotting.sessions.mlb(:,2);
tactileRod.mean = results.plotting.sessions.trb(:,2);
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));
    landmark.(sprintf('%s', session)) = results.plotting.sessions.lm(:,i+2);
    lineBisection.(sprintf('%s', session)) = results.plotting.sessions.mlb(:,i+2);
    tactileRod.(sprintf('%s', session)) = results.plotting.sessions.trb(:,i+2);
end

%% Individual tasks
cd(resDir)
% Number of bins - landmark
lmRange = max(landmark.mean) - min(landmark.mean);
nrBins = round(lmRange/0.5);
% standard deviation of population
lmSD = std(landmark.mean);
% Line bisection
mlbSD = std(lineBisection.mean);

%% Plotting landmark histo
pdfFileName = 'landmarkHisto_all.pdf';
figure()
lmH = histfit(landmark.mean, nrBins, 'kernel'); %histogram, 1 bin every 0.5mm, kernel distribution fit
% making pretty
lmH(1).FaceColor = [0.9 0.7 1]; 
lmH(2).Color = [0.2 0.1 0.3];
set(lmH(2), 'LineWidth', 3);
% axes
ax = gca;
lmHT = title('Landmark mean distribution'); set(lmHT, 'FontSize', 12);
lmHx = xlabel('Mean bias (mm)'); lmHy = ylabel('Frequency');
set(lmHx, 'FontSize', 12); set(lmHy, 'FontSize', 12);
set(ax, 'FontSize', 11, 'xtick', -6:1:6, 'ytick', 0:1:6);
yData = ax.YLim(1):ax.YLim(end); %adding 1 so can see end of histogram
xData = ax.XLim(1):ax.XLim(end);
shadedVal = zeros(1, length(yData));
% draw line at 0
hold on
lmHline = line('XData', [0,0], 'YData', [0 length(yData)], 'LineStyle', '-', ...
    'LineWidth', 1.5, 'Color', 'k');
% lines around 2*SD
hold on
lmHSD1 = line('XData', [lmSD*2 lmSD*2], 'YData', [0 length(yData)], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color', 'k');
hold on 
lmHSD2 = line('XData', [-lmSD*2 -lmSD*2], 'YData', [0 length(yData)], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color', 'k');
saveas(gcf, pdfFileName)

%% Plotting MLB histo
%% continue editing here - need an overall mean line too
pdfFileName = 'lineBisectionHisto_all.pdf';
figure()
mlbH = histfit(lineBisection.mean, nrBins, 'kernel'); %histogram, 1 bin every 0.5mm, kernel distribution fit
% making pretty
mlbH(1).FaceColor = [0.6 0.6 1]; 
mlbH(2).Color = [0.1 0.1 0.5];
set(mlbH(2), 'LineWidth', 3);
% axes
ax = gca;
mlbHT = title('Line bisection mean distribution'); set(mlbHT, 'FontSize', 12);
mlbHx = xlabel('Mean bias (mm)'); mlbHy = ylabel('Frequency');
set(mlbHx, 'FontSize', 12); set(mlbHy, 'FontSize', 12);
xlim([-12 12])
set(ax, 'FontSize', 11, 'xtick', [-12:2:12], 'ytick', 0:1:6);
yData = ax.YLim(1):ax.YLim(end); %adding 1 so can see end of histogram
xData = ax.XLim(1):ax.XLim(end);
shadedVal = zeros(1, length(yData));
% draw line at 0
hold on
mlbHline = line('XData', [0,0], 'YData', [0 length(yData)], 'LineStyle', '-', ...
    'LineWidth', 1.5, 'Color', 'k');
% lines around 2*SD
hold on
mlbHSD1 = line('XData', [mlbSD*2 mlbSD*2], 'YData', [0 length(yData)], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color', 'k');
hold on 
mlbHSD2 = line('XData', [-mlbSD*2 -mlbSD*2], 'YData', [0 length(yData)], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color', 'k');
saveas(gcf, pdfFileName)

%% Individual sessions - each task
%% Save and close