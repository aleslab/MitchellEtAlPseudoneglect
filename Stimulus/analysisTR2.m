%% Analysing Tactile Rod 2AFC task - AG.Mitchell 01.10.18
%% Data analysis script
% This script takes the data from the TR2 task and calcualtes the
% percentage of times participants said that the top line is longer
% This is then calculated against the amount (in mm) of shift in the top
% line (right or left)
% Psychometric curve calculated, and the point of subjective equality (PSE)
% used as a measure of bias in each participant 
% Psychometric function calculated using palamedes toolbox (type used:
% cumulative normal)

%% Note
% Lapse rate assumed 0 here, but will be calculated more accurately into
% paradigm 

%% Importing data
matfilename = ('TR2_results.mat');

data = struct;
data.size = [];
data.shift = [];
data.line = [];
data.response = [];
data.matrix = [];

[Data,Text] = xlsread('TrialMatrix_TR2'); %importing trial information for landmark tasks
data.size = Data(:, ismember(Text(1,:), 'size'));
data.shift = Data(:, ismember(Text(1,:), 'shift')); %shift in mm (0,2,4,6,8,1(=10))
data.line = Data(:, ismember(Text(1,:), 'line')); %indicating the line that should be shifted left (the other right, 1 = top, 2 = bottom)
data.response = Data(:, ismember(Text(1,:), 'response')); %top line = 1, bottom line = 2

%% Calculating percentage 'top line'
data.matrix = [data.size, data.shift, data.line, data.response];
top = data.matrix(find(data.matrix(:,4)==1),:); %all 'top longer' trials
bottom = data.matrix(find(data.matrix(:,4)==2),:); %all 'bottom longer' trials
topleft = data.matrix(find(data.matrix(:,4)==1 & data.matrix(:,3)==1),:);
bottomleft = data.matrix(find(data.matrix(:,4)==2 & data.matrix(:,3)==1),:); 
% calculating percentage of trials that the left shift line is thought to
% be longer
data.topl = (length(topleft)/length(top))*100;
data.bottoml = (length(bottomleft)/length(bottom))*100;
data.percentl = ((length(topleft)+length(bottomleft))/length(data.response))*100;

save(matfilename)
