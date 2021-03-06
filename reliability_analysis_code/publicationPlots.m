%% Bias experiment - poster plots
% Creating script that makes plots for A0 poster

%% Directories and data
% Getting directory
filePath = cd;
[dirBias, name, ext] = fileparts(filePath); %subject to change depending on where you analysis;
dirAna = [dirBias filesep 'Analysis']; %directory for all participant analysis
cd(dirAna)
load('ReliabilityAnalysis.mat');
mkdir Plots
dirPlot = [dirAna filesep 'Plots'];
cd(dirPlot)

%% Plotting - large scaled-up plots for posters
% creating figure for session subplots
pngFileName = strcat('Figure5', '.png');
figure('units', 'centimeters', 'Position', [10 1 17 55])

%colours for figures
s1_col = [0 0 0]; 
s2_col = [0.2 0.2 0.2];
s3_col = [0.4 0.4 0.4];
s4_col = [0.6 0.6 0.6];
mean_faceCol = [1 1 1];
mean_edgeCol = [0.1 0.1 0.1];

% Landmark task - to figure
% standard deviation values for shading
lmSDpt5 = allData.means.lmPSE(2)*0.5;
lmSD2 = allData.means.lmPSE(2)*2;
lmSDall = results.plotting.sessions.lm(:,7);
lmCI = results.plotting.sessions.lm(:,8); %need to do this again because data has been sorted

subplot(3,1,1)
% Adding SD shaded area - before the actual data for ease of visualisation
line('XData', [0 length(results.observers)], 'YData', [0, 0], 'LineStyle', '-', ...
    'LineWidth', 1, 'Color', 'k'); %midpoint
hold on
ax = gca;
xVal = [ax.XLim(1):ax.XLim(end)];
lm1 = scatter(results.plotting.sessions.lm(:,9), results.plotting.sessions.lm(:,3), ...
    'filled', 'o', 'MarkerFaceColor', s1_col); % session1
set(lm1, 'SizeData', 60);
hold on
lm2 = scatter(results.plotting.sessions.lm(:,9), results.plotting.sessions.lm(:,4), ...
    'filled', 'square', 'MarkerFaceColor', s2_col); % session1
set(lm2, 'SizeData', 60);
hold on
lm3 = scatter(results.plotting.sessions.lm(:,9), results.plotting.sessions.lm(:,5), ...
    'filled', 'o', 'MarkerFaceColor', s3_col); % session3
set(lm3, 'SizeData', 60);
hold on
lm4 = scatter(results.plotting.sessions.lm(:,9), results.plotting.sessions.lm(:,6), ...
    'filled', 'square', 'MarkerFaceColor', s4_col); % session4
set(lm4, 'SizeData', 60);
% Adding error bars to the mean data
hold on
lm5 = scatter(results.plotting.sessions.lm(:,9), results.plotting.sessions.lm(:,2), ...
    'filled', '^', 'MarkerFaceColor', mean_faceCol, 'MarkerEdgeColor', mean_edgeCol); % mean task data
set(lm5, 'SizeData', 80);
hold on
errorbar(results.plotting.sessions.lm(:,9), results.plotting.sessions.lm(:,2), lmCI, 'LineStyle', 'none',...
    'LineWidth', 1.5, 'Color', [0.1 0.1 0.1], 'CapSize', 0); %SD error bars first
ylim([-10 10]);
% Making it prettier
set(ax, 'FontSize', 11);
xLabels = num2str(results.plotting.sessions.lm(:,1));
xticks(ax, 1:length(results.plotting.sessions.lm(:,1)));
xticklabels(ax, xLabels);
ylabel('Bisection error (mm)');
%[lmlgd] = legend([lm1, lm2, lm3, lm4, lm5], '1', '2', '3', '4', 'Mean');
%lmText = [lmlgd, lmlgd.ItemText]; set(lmText, 'FontSize', 12);
%legend boxoff
box off
lmT = title('Landmark'); set(lmT, 'fontsize', 12);
%saveas(gcf, pngFileName);

% MLB - manual line bisection to figure
mlbSDpt5 = allData.means.mlb(2)*0.5;
mlbSD2 = allData.means.mlb(2)*2;
mlbSDall = results.plotting.sessions.mlb(:,7);
mlbCI = results.plotting.sessions.mlb(:,8); %need to do this again because data has been sorted

subplot(3,1,2)
% Adding SD shaded area - before the actual data for ease of visualisation
line('XData', [0 length(results.observers)], 'YData', [0, 0], 'LineStyle', '-', ...
    'LineWidth', 1, 'Color', 'k'); %midpoint
hold on
ax = gca;
xVal = [ax.XLim(1):ax.XLim(end)];
shadedVal = zeros(length(xVal),1); %making the same length so can plot shaded error bar around 0
mlb1 = scatter(results.plotting.sessions.mlb(:,9), results.plotting.sessions.mlb(:,3), ...
    'filled', 'o', 'MarkerFaceColor', s1_col); % session1
set(mlb1, 'SizeData', 60);
hold on
mlb2 = scatter(results.plotting.sessions.mlb(:,9), results.plotting.sessions.mlb(:,4), ...
    'filled', 'square', 'MarkerFaceColor', s2_col); % session1
set(mlb2, 'SizeData', 60);
hold on
mlb3 = scatter(results.plotting.sessions.mlb(:,9), results.plotting.sessions.mlb(:,5), ...
    'filled', 'o', 'MarkerFaceColor', s3_col); % session3
set(mlb3, 'SizeData', 60);
hold on
mlb4 = scatter(results.plotting.sessions.mlb(:,9), results.plotting.sessions.mlb(:,6), ...
    'filled', 'square', 'MarkerFaceColor', s4_col); % session4
set(mlb4, 'SizeData', 60);
% Adding error bars to the mean data
hold on
mlb5 = scatter(results.plotting.sessions.mlb(:,9), results.plotting.sessions.mlb(:,2), ...
    'filled', '^', 'MarkerFaceColor', mean_faceCol, 'MarkerEdgeColor', mean_edgeCol); % mean task data
set(mlb5, 'SizeData', 80);
errorbar(results.plotting.sessions.mlb(:,9), results.plotting.sessions.mlb(:,2), mlbCI, 'LineStyle', 'none',...
    'LineWidth', 1.5, 'Color', [0.1 0.1 0.1], 'CapSize', 0); %SD error bars first
hold on
ylim([-25 25]);
% Making it prettier
set(ax, 'FontSize', 11, 'ytick', [-20:10:20]);
xLabels = num2str(results.plotting.sessions.mlb(:,1));
xticks(ax, 1:length(results.plotting.sessions.mlb(:,1)));
xticklabels(ax, xLabels);ylabel('Bisection error (mm)');
box off
mlbT = title('Line bisection'); set(mlbT, 'fontsize', 12);
%saveas(gcf, pngFileName);

% TRB - adding to figure
% standard deviation values for shading
trbSDpt5 = allData.means.trb(2)*0.5;
trbSD2 = allData.means.trb(2)*2;
trbSDall = results.plotting.sessions.trb(:,7);
trbCI = results.plotting.sessions.trb(:,8); %need to do this again because data has been sorted

subplot(3,1,3)
% Adding SD shaded area - before the actual data for ease of visualisation

line('XData', [0 length(results.observers)], 'YData', [0, 0], 'LineStyle', '-', ...
    'LineWidth', 1, 'Color', 'k'); %midpoint
hold on
ax = gca;
xVal = [ax.XLim(1):ax.XLim(end)];
shadedVal = zeros(length(xVal),1); %making the same length so can plot shaded error bar around 0
trb1 = scatter(results.plotting.sessions.trb(:,9), results.plotting.sessions.trb(:,3), ...
    'filled', 'o', 'MarkerFaceColor', s1_col); % session1
set(trb1, 'SizeData', 60);
hold on
trb2 = scatter(results.plotting.sessions.trb(:,9), results.plotting.sessions.trb(:,4), ...
    'filled', 'square', 'MarkerFaceColor', s2_col); % session1
set(trb2, 'SizeData', 60);
hold on
trb3 = scatter(results.plotting.sessions.trb(:,9), results.plotting.sessions.trb(:,5), ...
    'filled', 'o', 'MarkerFaceColor', s3_col); % session3
set(trb3, 'SizeData', 60);
hold on
trb4 = scatter(results.plotting.sessions.trb(:,9), results.plotting.sessions.trb(:,6), ...
    'filled', 'square', 'MarkerFaceColor', s4_col); % session4
set(trb4, 'SizeData', 60);
% Adding error bars to the mean data
hold on
trb5 = scatter(results.plotting.sessions.trb(:,9), results.plotting.sessions.trb(:,2), ...
    'filled', '^', 'MarkerFaceColor', mean_faceCol, 'MarkerEdgeColor', mean_edgeCol); % mean task data
set(trb5, 'SizeData', 90);
errorbar(results.plotting.sessions.trb(:,9), results.plotting.sessions.trb(:,2), trbCI, 'LineStyle', 'none',...
    'LineWidth', 1.5, 'Color', [0.1 0.1 0.1], 'CapSize', 0); %SD error bars first
ylim([-25 25]);
% Making it prettier
set(ax, 'FontSize', 11, 'ytick', [-20:10:20]);
xLabels = num2str(results.plotting.sessions.trb(:,1));
xticks(ax, 1:length(results.plotting.sessions.trb(:,1)));
xticklabels(ax, xLabels);
xlabel('Observers'); ylabel('Bisection error (mm)');
%legend boxoff
% Adding text to define bias grouping
box off
trbT = title('Tactile rod bisection'); set(trbT, 'fontsize', 12);

% moving subplot legend
lgd = legend([trb1, trb2, trb3, trb4, trb5], '1', '2', '3', '4');
lgdText = [lgd, lgd.ItemText]; set(lgdText, 'FontSize', 11);
newLegendPos = [0.90 0.85 0.1 0.05];
newUnits = 'normalized';
set(lgd, 'Position', newLegendPos, 'Units', newUnits);
title(lgd, 'Session');
legend boxoff

print -r300 %set to resolution (let's try this)
saveas(gcf, pngFileName);


%% Modalities hypothesis
% standard deviation values for shading
SDpt5 = allData.means.tot(2)*0.5;
SD2 = allData.means.tot(2)*2;
SDall = results.plotting.modalities(:,6);
CIall = results.plotting.modalities(:,7);

pngFileName = strcat('Figure7', '.png');

figure('units', 'centimeters', 'Position', [5 3 18 10])
line('XData', [0 length(results.observers)], 'YData', [0, 0], 'LineStyle', '-', ...
    'LineWidth', 1, 'Color', 'k'); %midpoint
hold on
ax = gca;
xVal = [ax.XLim(1):ax.XLim(end)];
shadedVal = zeros(length(xVal),1); %making the same length so can plot shaded error bar around 0
m1 = scatter(results.plotting.modalities(:,8), results.plotting.modalities(:,3), ...
    'filled', 'o', 'MarkerFaceColor', [0.3 0.3 0.3]); % landmark task data
set(m1, 'SizeData', 60);
hold on
m2 = scatter(results.plotting.modalities(:,8), results.plotting.modalities(:,4), ...
    'filled', 'o', 'MarkerFaceColor', [0.1 0.2 0.6]); % mlb task data
set(m2, 'SizeData', 60);
hold on
m3 = scatter(results.plotting.modalities(:,8), results.plotting.modalities(:,5), ...
    'filled', 'o', 'MarkerFaceColor', [0.1 0.5 0.2]); % trb task data
set(m3, 'SizeData', 60);
hold on
m4 = scatter(results.plotting.modalities(:,8), results.plotting.modalities(:,2), ...
    'filled', '^', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0.1 0.1 0.1]); % mean task data
hold on
errorbar(results.plotting.modalities(:,8), results.plotting.modalities(:,2), CIall, 'LineStyle', 'none',...
    'LineWidth', 1.5, 'Color', [0 0 0], 'CapSize', 0);
set(m4, 'SizeData', 80);
ylim([-15 15]);
% Adding SD shaded area
% Adding error bars to the mean data
% Making it prettier
set(ax, 'FontSize', 11);
xLabels = num2str(results.plotting.modalities(:,1));
xticks(ax, 1:length(results.plotting.modalities(:,1)));
xticklabels(ax, xLabels);
%set(ax, 'XTick', results.observers);
xlabel('Observers'); ylabel('Bisection error (mm)');
mlgd = legend([m1 m2 m3 m4], 'Landmarks', 'MLB', 'TRB', 'Mean', [110 230 0.2 0.1]);
legend boxoff
mText = [mlgd, mlgd.ItemText]; set(mText, 'FontSize', 10);
% Adding text to define bias grouping
%mT = title('Modalities'); set(mT, 'fontsize', 34);
print -r300
saveas(gcf, pngFileName);

%% Modalities plot showing varaition within participants
% Like the within data smooth pursuit plots
% Variables needed for the plot
nData = length(results.observers); 
barwidth = 1;
xJitter = 0.01*randn(nData,1);
xLocMatrix = [xJitter+.6 xJitter+2 xJitter+3.4]'; %Make a matrix of where to plot
dataMat = [results.plotting.modalities(:,3), results.plotting.modalities(:,4), ...
    results.plotting.modalities(:,5)]';

% Making plot
pngFileName = 'Figure8.png';
figure('units', 'centimeters', 'Position', [5 3 12 8]);
clf;set(gcf,'color','w');hold on;

%Draw invidivual datapoints and the lines
lineH= plot(xLocMatrix, dataMat,'-');
hold on
symbolH = plot(xLocMatrix, dataMat,'.'); 
%change the line to be light gray. 
set(lineH,'color',[.75 .75 .75], 'linewidth', 1.5)
%Set the points to be a bit darker and a nice size. 
set(symbolH,'color',[.3 .3 .3],'markersize',12);
%Just bring up a bar plot, make it
barH = bar([.5 2 3.5],[mean(dataMat(1,:)), mean(dataMat(2,:)), mean(dataMat(3,:))], barwidth)
axH = gca; %Get the current axes the bar got plotteed into. 
%Let's changeup the bar now.
set(barH,'facecolor','none','linestyle','-','linewidth',2, 'barwidth', 0.7);
%Next put some std error bars on the bar plot:
errH = errorbar([.5 2 3.5],[mean(dataMat(1,:)) mean(dataMat(2,:)) mean(dataMat(3,:))],...
    [std(dataMat(1,:))/sqrt(nData) std(dataMat(2,:))/sqrt(nData) std(dataMat(3,:))/sqrt(nData)],'.')
set(errH,'color','k','linewidth',1.5, 'capsize', 7);
%Setup the axis/background
set(axH,'linewidth',1,'fontsize',11,'xtick',[.5 2 3.5], 'ytick', [-8:2:6]);
%axH.TickLength = [.03 0.1];
axH.XTickLabel = {'landmark' 'line bisec.' 'rod bisec.'}; axH.YLim = [-9 7];
axH = ylabel(['Bisection error (mm)']);
%mT = title('Means all tasks'); set(mT, 'fontsize', 32);
%tH = title('Bias within modalities'); set(tH,'fontsize',30);
box off
print -r300
saveas(gcf, pngFileName);

%% Sessions plot showing variation within participants
nData = length(results.observers); 
barwidth = 1;
xJitter = 0.01*randn(nData,1);
xLocMatrix = [xJitter+.6 xJitter+1.6 xJitter+2.4 xJitter+3.4]'; %Make a matrix of where to plot
dataMat = allData.sessions.allSessions';

% Making plot
pngFileName = 'Figure6.png';
figure('units', 'centimeters', 'Position', [5 3 12 10]);
clf;set(gcf,'color','w');hold on;

%Draw invidivual datapoints and the lines
lineH= plot(xLocMatrix, dataMat,'-');
hold on
symbolH = plot(xLocMatrix, dataMat,'.'); 
%change the line to be light gray. 
set(lineH,'color',[.75 .75 .75], 'linewidth', 1.5)
%Set the points to be a bit darker and a nice size. 
set(symbolH,'color',[.3 .3 .3],'markersize',12);
%Just bring up a bar plot, make it
barH = bar([.5 1.5 2.5 3.5],[mean(dataMat(1,:)), mean(dataMat(2,:)), mean(dataMat(3,:)),...
    mean(dataMat(4,:))], barwidth);
axH = gca; %Get the current axes the bar got plotteed into. 
%Let's changeup the bar now.
set(barH,'facecolor','none','linestyle','-','linewidth',2, 'barwidth', 0.7);
%Next put some std error bars on the bar plot:
errH = errorbar([.5 1.5 2.5 3.5],[mean(dataMat(1,:)) mean(dataMat(2,:)) mean(dataMat(3,:)), mean(dataMat(4,:))],...
    [std(dataMat(1,:))/sqrt(nData) std(dataMat(2,:))/sqrt(nData)...
    std(dataMat(3,:))/sqrt(nData), std(dataMat(4,:))/sqrt(nData)],'.');
set(errH,'color','k','linewidth',1.5, 'capsize', 7);
%Setup the axis/background
set(axH,'linewidth',1,'fontsize',11,'xtick',[.5 1.5 2.5 3.5], 'ytick', [-8:2:6]);
%axH.TickLength = [.03 0.1];
axH.XTickLabel = {'session 1' ' sesion 2' 'session 3' 'session 4'}; axH.YLim = [-9 7];
axH = ylabel(['Bisection error (mm)']);
%mT = title('Means all tasks'); set(mT, 'fontsize', 32);
%tH = title('Bias within modalities'); set(tH,'fontsize',30);
box off
print -r300
saveas(gcf, pngFileName);



