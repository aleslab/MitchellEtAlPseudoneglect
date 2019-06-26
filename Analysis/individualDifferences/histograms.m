%% AG. Mitchell - 14.06.19
%% Plotting individual participant data into a histogram
% Do different tasks show different spread of data?
% Look acros individual tasks and also individual testing sessions

% Take data from the reliabilityAnalysis.mat file (extracted using
% ReliabilityAnalysis.m script)

%% Data
% pathways
expDir = 'M:\Alex_Files\Experiments\Bias\Analysis\';
dataDir = [expDir]; %original data set
resDir = [expDir 'individual_differences']; %for saving plots/individual diff analysis

% getting data
cd(dataDir)
load('reliabilityAnalysis.mat');
matFileName = 'distributionAnalysis.mat';

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
lmMean = mean(landmark.mean);
% Line bisection
mlbSD = std(lineBisection.mean);
mlbMean = mean(lineBisection.mean);
% Rod bisection
trbSD = std(tactileRod.mean);
trbMean = mean(tactileRod.mean);

%% Plotting landmark histo
pngFileName = 'histo_all.png';
figure('pos',[150 150 1600 400])
subplot(1,3,1)
% setting characteristics
lmH = histfit(landmark.mean, nrBins, 'kernel'); %histogram, 1 bin every 0.5mm, kernel distribution fit
% making pretty
lmH(1).FaceColor = [0.9 0.7 1]; 
lmH(2).Color = [0.2 0.1 0.3];
set(lmH(2), 'LineWidth', 3);
% axes
ax = gca;
set(ax, 'FontSize', 11, 'xtick', -6:1:6, 'ytick', 0:1:6);
yData = ax.YLim(1):ax.YLim(end); %adding 1 so can see end of histogram
xData = ax.XLim(1):ax.XLim(end);
lmHT = title('Landmark mean distribution'); set(lmHT, 'FontSize', 12);
lmHx = xlabel('Mean bias (mm)'); lmHy = ylabel('Frequency');
set(lmHx, 'FontSize', 12); set(lmHy, 'FontSize', 12);
shadedVal = zeros(1, length(yData));
% draw line at 0
hold on
lmHline = line('XData', [0,0], 'YData', [0 length(yData)], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color', 'k');
% draw line at mean
lmHmean = line('XData', [lmMean,lmMean], 'YData', [0 length(yData)], 'LineStyle', '-', ...
    'LineWidth', 2, 'Color', [0.2 0.1 0.3]);
% lines around 2*SD
lmHSD1 = line('XData', [lmSD*2 lmSD*2], 'YData', [0 length(yData)], 'LineStyle', '-', ...
    'LineWidth', 1, 'Color', [0.4 0.4 0.4]);
hold on 
lmHSD2 = line('XData', [-lmSD*2 -lmSD*2], 'YData', [0 length(yData)], 'LineStyle', '-', ...
    'LineWidth', 1, 'Color', [0.4 0.4 0.4]);
%saveas(gcf, pdfFileName)

%% Plotting MLB histo
%pdfFileName = 'lineBisectionHisto_all.pdf';
subplot(1,3,2)
mlbH = histfit(lineBisection.mean, nrBins, 'kernel'); %histogram, 1 bin every 0.5mm, kernel distribution fit
% making pretty
mlbH(1).FaceColor = [0.7 0.7 1]; 
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
mlbHline = line('XData', [0,0], 'YData', [0 length(yData)], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color', 'k');
% line at mean
mlbHmean = line('XData', [mlbMean,mlbMean], 'YData', [0 length(yData)], 'LineStyle', '-', ...
    'LineWidth', 2, 'Color', [0.1 0.1 0.5]);
% lines around 2*SD
hold on
mlbHSD1 = line('XData', [mlbSD*2 mlbSD*2], 'YData', [0 length(yData)], 'LineStyle', '-', ...
    'LineWidth', 1, 'Color', [0.4 0.4 0.4]);
hold on 
mlbHSD2 = line('XData', [-mlbSD*2 -mlbSD*2], 'YData', [0 length(yData)], 'LineStyle', '-', ...
    'LineWidth', 1, 'Color', [0.4 0.4 0.4]);
%saveas(gcf, pdfFileName)

%% Plotting TRB histo
subplot(1,3,3)
trbH = histfit(tactileRod.mean, nrBins, 'kernel'); %histogram, 1 bin every 0.5mm, kernel distribution fit
% making pretty
trbH(1).FaceColor = [1 0.7 0.8]; 
trbH(2).Color = [0.4 0 0.2];
set(trbH(2), 'LineWidth', 3);
% axes
ax = gca;
trbHT = title('Tactile rod mean distribution'); set(trbHT, 'FontSize', 12);
trbHx = xlabel('Mean bias (mm)'); trbHy = ylabel('Frequency');
set(trbHx, 'FontSize', 12); set(trbHy, 'FontSize', 12);
xlim([-12 12])
set(ax, 'FontSize', 11, 'xtick', [-12:2:12], 'ytick', 0:1:6);
yData = ax.YLim(1):ax.YLim(end); %adding 1 so can see end of histogram
xData = ax.XLim(1):ax.XLim(end);
shadedVal = zeros(1, length(yData));
% draw line at 0
hold on
trbHline = line('XData', [0,0], 'YData', [0 length(yData)], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color', 'k');
% line at mean
trbHmean = line('XData', [trbMean,trbMean], 'YData', [0 length(yData)], 'LineStyle', '-', ...
    'LineWidth', 2, 'Color', [0.4 0 0.2]);
% lines around 2*SD
hold on
trbHSD1 = line('XData', [trbSD*2 trbSD*2], 'YData', [0 length(yData)], 'LineStyle', '-', ...
    'LineWidth', 1, 'Color', [0.4 0.4 0.4]);
hold on 
trbHSD2 = line('XData', [-trbSD*2 -trbSD*2], 'YData', [0 length(yData)], 'LineStyle', '-', ...
    'LineWidth', 1, 'Color', [0.4 0.4 0.4]);
% h = gcf; set(h, 'PaperOrientation', 'landscape');
% set(h, 'PaperPositionMode', 'auto'); set(h, 'PaperPosition', [1 1 30 14]);
saveas(gcf, pngFileName)

%% K-S testing to check for normality
% one sample ks test to check the spread of each data-set
[normStat.lm.h, normStat.lm.p] = kstest(landmark.mean); %landmark
[normStat.mlb.h, normStat.mlb.p] = kstest(lineBisection.mean); %line bisection
[normStat.trb.h, normStat.trb.p] = kstest(tactileRod.mean); %tactile rod 

%% Bayesian gaussian fit test
% Code by J Ales - 24.06.19
% Landmark task gaussian fit
X = landmark.mean; %Same data you put into the histogram.  
gmOnePop = fitgmdist(X,1); %Fit 1 gaussian to distribution
gmTwoPop = fitgmdist(X,2); %Fit mixture of 2 gaussians to distribution
gmFit.landmark.onePop = gmOnePop; %adding info to own struct 
gmFit.landmark.twoPop = gmTwoPop;
%Compare the evidence:
gmFit.landmark.bayesFactor = gmTwoPop.BIC - gmOnePop.BIC; % (+) means in favor of ONE, (-) means in favor of Two

% Line bisection gaussian fit
X = lineBisection.mean; %Same data you put into the histogram.  
gmOnePop = fitgmdist(X,1); %Fit 1 gaussian to distribution
gmTwoPop = fitgmdist(X,2); %Fit mixture of 2 gaussians to distribution
gmFit.lineBisection.onePop = gmOnePop; %adding info to own struct 
gmFit.lineBisection.twoPop = gmTwoPop;
%Compare the evidence:
gmFit.lineBisection.bayesFactor = gmTwoPop.BIC - gmOnePop.BIC; % (+) means in favor of ONE, (-) means in favor of Two

% Tactile rod gaussian fit
X = tactileRod.mean; %Same data you put into the histogram.  
gmOnePop = fitgmdist(X,1); %Fit 1 gaussian to distribution
gmTwoPop = fitgmdist(X,2); %Fit mixture of 2 gaussians to distribution
gmFit.tactileRod.onePop = gmOnePop; %adding info to own struct 
gmFit.tactileRod.twoPop = gmTwoPop;
%Compare the evidence:
gmFit.tactileRod.bayesFactor = gmTwoPop.BIC - gmOnePop.BIC; % (+) means in favor of ONE, (-) means in favor of Two

%% Individual sessions - each task
%% Save and close
save(matFileName, 'gmFit', 'normStat', 'landmark', 'lineBisection', 'tactileRod');