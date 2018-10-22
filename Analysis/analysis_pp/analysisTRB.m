%% Analysis for Tactile Rod Bisection - AG.Mitchell 01.10.18
%% Runs analysis on TRB task
% Loads in the data for the task and averages the response across all
% variables
% Then calculates bias
% Then calculates bias relative to line length and shift into L/R hemifield

%% Importing the data
ppID = input('Participant ID? ', 's'); %for use when navigating files
matfilename = sprintf('TRB_analysis%s.mat', ppID);

data = struct;
data.size = [];
data.side = [];
data.response = [];
data.matrix = [];

[Data,Text] = xlsread('TrialMatrix_TRB'); %importing trial information for landmark tasks
data.size = Data(:, ismember(Text(1,:), 'size'));
data.side = Data(:, ismember(Text(1,:), 'shift')); %2deg left = 1, middle = 0, 2deg right = 2
data.response = Data(:, ismember(Text(1,:), 'response')); %in mm, relative to length of line

%% Finding error 
data.size = data.size*100; %finding actual line length in mm to calculate error
data.error = data.response - data.size/2; %calculating data error in mm
% data matrix for later analyses
data.matrix = [data.size, data.side, data.response, data.error];
data.errormean = nanmean(data.error); %average results in mm
data.errorsd = nanstd(data.error);

save(matfilename);