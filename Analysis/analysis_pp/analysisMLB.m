%% A.G. Mitchell 22.10.18
%% Imports MLB data from each session for analysis
% This is the analysis script run on each participant
% Takes the average error from each session and calculates total
% participant bias
% Also important to maintain mean biases in each session for comparison in
% the second (across time) hypothesis
% Places data in a format readable for the overall analysis

%% Variables
ppID = input('Participant ID? ', 's'); %for use when navigating files
matfilename = sprintf('MLB_analysis%s.mat', ppID);

% Directory
dir_data = 'C:\Users\Experimenter\Documents\Experiments2018\Bias\Data';
dir_PP = [dir_data filesep ppID filesep];


%% Taking data from each session
%% Average values