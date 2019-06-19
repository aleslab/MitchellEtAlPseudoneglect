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

%% Landmark vs. LM2AFC
