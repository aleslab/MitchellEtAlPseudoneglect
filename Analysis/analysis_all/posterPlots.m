%% Bias experiment - poster plots
% Creating script that makes plots for A0 poster

%% Directories and data
dirBias = ('M:\Experiments\Bias');
dirAna = [dirBias filesep 'Analysis']; %directory for all participant analysis
cd(dirAna)
load('ReliabilityAnalysis.mat');
plotDir = ('M:\Presentations\VSS2019\Figures');
cd(plotDir)

%% Plotting - large scaled-up plots for posters
%% Landmark task
% standard deviation values for shading
lmSDpt5 = allData.means.lmPSE(2)*0.5;
lmSD2 = allData.means.lmPSE(2)*2;
lmSDall = results.plotting.sessions.lm(:,7);
lmCI = results.plotting.sessions.lm(:,8); %need to do this again because data has been sorted

pngFileName = strcat('biasSessions_lm', '.png');

figure('units', 'centimeters', 'Position', [5 3 34 18])
errorbar(results.plotting.sessions.lm(:,9), results.plotting.sessions.lm(:,2), lmSDall, 'LineStyle', 'none',...
    'LineWidth', 1.5, 'Color', [0.1 0.1 0.1], 'CapSize', 0); %SD error bars first
hold on
% Adding SD shaded area - before the actual data for ease of visualisation
ax = gca;
xVal = [ax.XLim(1):ax.XLim(end)];
shadedVal = zeros(length(xVal),1); %making the same length so can plot shaded error bar around 0
createShadedRegion(xVal, shadedVal, (shadedVal - lmSDpt5), (shadedVal + lmSDpt5),':','color', [0.3 0.3 0.3]);
hold on
createShadedRegion(xVal, shadedVal, (shadedVal - lmSD2), (shadedVal + lmSD2),':','color', [0.5 0.5 0.5]);
hold on
line('XData', [0 length(results.observers)], 'YData', [0, 0], 'LineStyle', '-', ...
    'LineWidth', 0.5, 'Color', 'k'); %midpoint
hold on
lm1 = scatter(results.plotting.sessions.lm(:,9), results.plotting.sessions.lm(:,3), ...
    'filled', 'o', 'MarkerFaceColor', [0.1 0 0.2]); % session1
set(lm1, 'SizeData', 90);
hold on
lm2 = scatter(results.plotting.sessions.lm(:,9), results.plotting.sessions.lm(:,4), ...
    'filled', 'o', 'MarkerFaceColor', [0.3 0.1 0.4]); % session1
set(lm2, 'SizeData', 90);
hold on
lm3 = scatter(results.plotting.sessions.lm(:,9), results.plotting.sessions.lm(:,5), ...
    'filled', 'o', 'MarkerFaceColor', [0.5 0.2 0.6]); % session3
set(lm3, 'SizeData', 90);
hold on
lm4 = scatter(results.plotting.sessions.lm(:,9), results.plotting.sessions.lm(:,6), ...
    'filled', 'o', 'MarkerFaceColor', [0.7 0.5 0.8]); % session4
set(lm4, 'SizeData', 90);
% Adding error bars to the mean data
hold on
lm5 = scatter(results.plotting.sessions.lm(:,9), results.plotting.sessions.lm(:,2), ...
    'filled', '^', 'MarkerFaceColor', [0.9 0.8 1], 'MarkerEdgeColor', [0.3 0 0.4]); % mean task data
set(lm5, 'SizeData', 100);
ylim([-10 10]);
% Making it prettier
set(ax, 'FontSize', 20);
xLabels = num2str(results.plotting.sessions.lm(:,1));
xticks(ax, 1:length(results.plotting.sessions.lm(:,1)));
xticklabels(ax, xLabels);
xlabel('Observers'); ylabel('Bias (mm)');
lmlgd = legend([lm1, lm2, lm3, lm4, lm5], '1', '2', '3', '4', 'Mean', [175 400 0.1 0.2]);
lmText = [lmlgd, lmlgd.ItemText]; set(lmText, 'FontSize', 20);
legend boxoff
% Adding text to define bias grouping
leftDim = [0.17 0.13 0.235 0.045];
rightDim = [0.575 0.13 0.29 0.045]; midDim = [0.41 0.13 0.16 0.045];
% annotation('textbox', leftDim, 'String', 'Left', 'FontSize', 8); %'Color', [0.2 0.2 0.2]);
% annotation('textbox', rightDim, 'String', 'Right', 'FontSize', 8); %'Color', [0.7 0.7 0.7]);
% annotation('textbox', midDim, 'String', 'Unspecified', 'FontSize', 8);  %'Color', [0.45 0.45 0.45]);
% % saving as both pdf and png for ease
box off
lmT = title('Landmark'); set(lmT, 'fontsize', 30);
saveas(gcf, pngFileName);

%% MLB
mlbSDpt5 = allData.means.mlb(2)*0.5;
mlbSD2 = allData.means.mlb(2)*2;
mlbSDall = results.plotting.sessions.mlb(:,7);
mlbCI = results.plotting.sessions.mlb(:,8); %need to do this again because data has been sorted

pngFileName = strcat('biasSessions_mlb', '.png');

figure('units', 'centimeters', 'Position', [5 3 34 18])
errorbar(results.plotting.sessions.mlb(:,9), results.plotting.sessions.mlb(:,2), mlbSDall, 'LineStyle', 'none',...
    'LineWidth', 1.5, 'Color', [0.1 0.1 0.1], 'CapSize', 0); %SD error bars first
hold on
% Adding SD shaded area - before the actual data for ease of visualisation
ax = gca;
xVal = [ax.XLim(1):ax.XLim(end)];
shadedVal = zeros(length(xVal),1); %making the same length so can plot shaded error bar around 0
createShadedRegion(xVal, shadedVal, (shadedVal - mlbSDpt5), (shadedVal + mlbSDpt5),':','color', [0.3 0.3 0.3]);
hold on
createShadedRegion(xVal, shadedVal, (shadedVal - mlbSD2), (shadedVal + mlbSD2),':','color', [0.5 0.5 0.5]);
hold on
line('XData', [0 length(results.observers)], 'YData', [0, 0], 'LineStyle', '-', ...
    'LineWidth', 0.5, 'Color', 'k'); %midpoint
hold on
mlb1 = scatter(results.plotting.sessions.mlb(:,9), results.plotting.sessions.mlb(:,3), ...
    'filled', 'o', 'MarkerFaceColor', [0 0.1 0.2]); % session1
set(mlb1, 'SizeData', 90);
hold on
mlb2 = scatter(results.plotting.sessions.mlb(:,9), results.plotting.sessions.mlb(:,4), ...
    'filled', 'o', 'MarkerFaceColor', [0.1 0.2 0.5]); % session1
set(mlb2, 'SizeData', 90);
hold on
mlb3 = scatter(results.plotting.sessions.mlb(:,9), results.plotting.sessions.mlb(:,5), ...
    'filled', 'o', 'MarkerFaceColor', [0.1 0.4 0.7]); % session3
set(mlb3, 'SizeData', 90);
hold on
mlb4 = scatter(results.plotting.sessions.mlb(:,9), results.plotting.sessions.mlb(:,6), ...
    'filled', 'o', 'MarkerFaceColor', [0.4 0.6 0.9]); % session4
set(mlb4, 'SizeData', 90);
% Adding error bars to the mean data
hold on
mlb5 = scatter(results.plotting.sessions.mlb(:,9), results.plotting.sessions.mlb(:,2), ...
    'filled', '^', 'MarkerFaceColor', [0.6 0.8 1], 'MarkerEdgeColor', [0 0.1 0.3]); % mean task data
set(mlb5, 'SizeData', 100);
ylim([-25 25]);
% Making it prettier
set(ax, 'FontSize', 20, 'ytick', [-20:10:20]);
xLabels = num2str(results.plotting.sessions.mlb(:,1));
xticks(ax, 1:length(results.plotting.sessions.mlb(:,1)));
xticklabels(ax, xLabels);
xlabel('Observers'); ylabel('Bias (mm)');
mlblgd = legend([mlb1, mlb2, mlb3, mlb4, mlb5], '1', '2', '3', '4', 'Mean', [175 400 0.1 0.2]);
mlbText = [mlblgd, mlblgd.ItemText]; set(mlbText, 'FontSize', 20);
legend boxoff
% Adding text to define bias grouping
leftDim = [0.17 0.13 0.235 0.045];
rightDim = [0.575 0.13 0.29 0.045]; midDim = [0.41 0.13 0.16 0.045];
% annotation('textbox', leftDim, 'String', 'Left', 'FontSize', 8); %'Color', [0.2 0.2 0.2]);
% annotation('textbox', rightDim, 'String', 'Right', 'FontSize', 8); %'Color', [0.7 0.7 0.7]);
% annotation('textbox', midDim, 'String', 'Unspecified', 'FontSize', 8);  %'Color', [0.45 0.45 0.45]);
% % saving as both pdf and png for ease
box off
mlbT = title('Line bisection'); set(mlbT, 'fontsize', 30);
saveas(gcf, pngFileName);

%% TRB
% standard deviation values for shading
trbSDpt5 = allData.means.trb(2)*0.5;
trbSD2 = allData.means.trb(2)*2;
trbSDall = results.plotting.sessions.trb(:,7);
trbCI = results.plotting.sessions.trb(:,8); %need to do this again because data has been sorted

pngFileName = strcat('biasSessions_trb', '.png');

figure('units', 'centimeters', 'Position', [5 3 34 18])
errorbar(results.plotting.sessions.trb(:,9), results.plotting.sessions.trb(:,2), trbSDall, 'LineStyle', 'none',...
    'LineWidth', 1.5, 'Color', [0.1 0.1 0.1], 'CapSize', 0); %SD error bars first
hold on
% Adding SD shaded area - before the actual data for ease of visualisation
ax = gca;
xVal = [ax.XLim(1):ax.XLim(end)];
shadedVal = zeros(length(xVal),1); %making the same length so can plot shaded error bar around 0
createShadedRegion(xVal, shadedVal, (shadedVal - trbSDpt5), (shadedVal + trbSDpt5),':','color', [0.3 0.3 0.3]);
hold on
createShadedRegion(xVal, shadedVal, (shadedVal - trbSD2), (shadedVal + trbSD2),':','color', [0.5 0.5 0.5]);
hold on
line('XData', [0 length(results.observers)], 'YData', [0, 0], 'LineStyle', '-', ...
    'LineWidth', 0.5, 'Color', 'k'); %midpoint
hold on
trb1 = scatter(results.plotting.sessions.trb(:,9), results.plotting.sessions.trb(:,3), ...
    'filled', 'o', 'MarkerFaceColor', [0.1 0.1 0.2]); % session1
set(trb1, 'SizeData', 90);
hold on
trb2 = scatter(results.plotting.sessions.trb(:,9), results.plotting.sessions.trb(:,4), ...
    'filled', 'o', 'MarkerFaceColor', [0.2 0.2 0.4]); % session1
set(trb2, 'SizeData', 90);
hold on
trb3 = scatter(results.plotting.sessions.trb(:,9), results.plotting.sessions.trb(:,5), ...
    'filled', 'o', 'MarkerFaceColor', [0.5 0.5 0.7]); % session3
set(trb3, 'SizeData', 90);
hold on
trb4 = scatter(results.plotting.sessions.trb(:,9), results.plotting.sessions.trb(:,6), ...
    'filled', 'o', 'MarkerFaceColor', [0.7 0.7 0.9]); % session4
set(trb4, 'SizeData', 90);
% Adding error bars to the mean data
hold on
trb5 = scatter(results.plotting.sessions.trb(:,9), results.plotting.sessions.trb(:,2), ...
    'filled', '^', 'MarkerFaceColor', [0.9 0.9 1], 'MarkerEdgeColor', [0.1 0.1 0.2]); % mean task data
set(trb5, 'SizeData', 100);
ylim([-25 25]);
% Making it prettier
set(ax, 'FontSize', 20, 'ytick', [-20:10:20]);
xLabels = num2str(results.plotting.sessions.trb(:,1));
xticks(ax, 1:length(results.plotting.sessions.trb(:,1)));
xticklabels(ax, xLabels);
xlabel('Observers'); ylabel('Bias (mm)');
trblgd = legend([trb1, trb2, trb3, trb4, trb5], '1', '2', '3', '4', 'Mean', [175 400 0.1 0.2]);
trbText = [trblgd, trblgd.ItemText]; set(trbText, 'FontSize', 20);
legend boxoff
% Adding text to define bias grouping
leftDim = [0.17 0.13 0.235 0.045];
rightDim = [0.575 0.13 0.29 0.045]; midDim = [0.41 0.13 0.16 0.045];
% annotation('textbox', leftDim, 'String', 'Left', 'FontSize', 8); %'Color', [0.2 0.2 0.2]);
% annotation('textbox', rightDim, 'String', 'Right', 'FontSize', 8); %'Color', [0.7 0.7 0.7]);
% annotation('textbox', midDim, 'String', 'Unspecified', 'FontSize', 8);  %'Color', [0.45 0.45 0.45]);
% % saving as both pdf and png for ease
box off
trbT = title('Tactile rod bisection'); set(trbT, 'fontsize', 30);
saveas(gcf, pngFileName);

%% All sessions
% standard deviation values for shading
allSDpt5 = std(allData.means.allPP(:,1))*0.5;
allSD2 = std(allData.means.allPP(:,1))*2;
allSDall = results.plotting.sessions.all(:,7);
allCI = results.plotting.sessions.all(:,8); %need to do this again because data has been sorted
pngFileName = strcat('biasSessions_all', '.png');

figure('units', 'centimeters', 'Position', [5 3 34 18])
errorbar(results.plotting.sessions.all(:,9), results.plotting.sessions.all(:,2), allSDall, 'LineStyle', 'none',...
    'LineWidth', 1.5, 'Color', [0.1 0.1 0.1], 'CapSize', 0); %SD error bars first
hold on
% Adding SD shaded area - before the actual data for ease of visualisation
ax = gca;
xVal = [ax.XLim(1):ax.XLim(end)];
shadedVal = zeros(length(xVal),1); %making the same length so can plot shaded error bar around 0
createShadedRegion(xVal, shadedVal, (shadedVal - allSDpt5), (shadedVal + allSDpt5),':','color', [0.3 0.3 0.3]);
hold on
createShadedRegion(xVal, shadedVal, (shadedVal - allSD2), (shadedVal + allSD2),':','color', [0.7 0.7 0.7]);
hold on
line('XData', [0 length(results.observers)], 'YData', [0, 0], 'LineStyle', '-', ...
    'LineWidth', 0.5, 'Color', 'k'); %midpoint
hold on
all1 = scatter(results.plotting.sessions.all(:,9), results.plotting.sessions.all(:,3), ...
    'filled', 'o', 'MarkerFaceColor', [0 0 0]); % session1
set(all1, 'SizeData', 90);
hold on
all2 = scatter(results.plotting.sessions.all(:,9), results.plotting.sessions.all(:,4), ...
    'filled', 'o', 'MarkerFaceColor', [0.2 0.2 0.2]); % session1
set(all2, 'SizeData', 90);
hold on
all3 = scatter(results.plotting.sessions.all(:,9), results.plotting.sessions.all(:,5), ...
    'filled', 'o', 'MarkerFaceColor', [0.4 0.4 0.4]); % session3
set(all3, 'SizeData', 90);
hold on
all4 = scatter(results.plotting.sessions.all(:,9), results.plotting.sessions.all(:,6), ...
    'filled', 'o', 'MarkerFaceColor', [0.6 0.6 0.6]); % session4
set(all4, 'SizeData', 90);
% Adding error bars to the mean data
hold on
all5 = scatter(results.plotting.sessions.all(:,9), results.plotting.sessions.all(:,2), ...
    'filled', '^', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerEdgeColor', [0.1 0.1 0.1]); % mean task data
set(all5, 'SizeData', 100);
ylim([-15 15]);
% Making it prettier
set(ax, 'FontSize', 20, 'ytick', [-15:5:15]);
xLabels = num2str(results.plotting.sessions.all(:,1));
xticks(ax, 1:length(results.plotting.sessions.all(:,1)));
xticklabels(ax, xLabels);
xlabel('Observers'); ylabel('Bias (mm)');
alllgd = legend([all1, all2, all3, all4, all5], '1', '2', '3', '4', 'Mean', [175 400 0.1 0.2]);
allText = [alllgd, alllgd.ItemText]; set(allText, 'FontSize', 20);
legend boxoff
% Adding text to define bias grouping
leftDim = [0.17 0.13 0.235 0.045];
rightDim = [0.575 0.13 0.29 0.045]; midDim = [0.41 0.13 0.16 0.045];
% annotation('textbox', leftDim, 'String', 'Left', 'FontSize', 8); %'Color', [0.2 0.2 0.2]);
% annotation('textbox', rightDim, 'String', 'Right', 'FontSize', 8); %'Color', [0.7 0.7 0.7]);
% annotation('textbox', midDim, 'String', 'Unspecified', 'FontSize', 8);  %'Color', [0.45 0.45 0.45]);
% % saving as both pdf and png for ease
box off
allT = title('All tasks'); set(allT, 'fontsize', 30);
saveas(gcf, pngFileName);

