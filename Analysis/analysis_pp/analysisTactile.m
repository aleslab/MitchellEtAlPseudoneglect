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

matfilename = sprintf('TR2_analysis%s.mat', ppID);

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
for i = 1:length(data.response)
    perceivedlong(i) = isequal(data.line(i), data.response(i));
end
data.percentl = (sum(perceivedlong)/length(data.response))*100;

save(matfilename)
