%% A.G. Mitchell 22.10.18
%% Imports MLB data from each session for analysis
% This is the analysis script run on each participant
% Takes the average error from each session and calculates total
% participant bias
% Also important to maintain mean biases in each session for comparison in
% the second (across time) hypothesis
% Places data in a format readable for the overall analysis

clear all

%% Variables
ppID = input('Participant ID? ', 's'); %for use when navigating files
matfilename = sprintf('MLB_analysis%s.mat', ppID);
nSessions = 1:4; %vector number of sessions each participant does
% Directory
dirBias = ('C:\Users\Experimenter\Documents\Experiments2018\Bias'); %subject to change depending on where you analyse
dirPP = [dirBias filesep ppID filesep]; %participant directory

% Load data for all trials
cd(dirPP)
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));
    dirSess = [dirPP filesep session filesep]; %current session pathway
    cd(dirSess); %directing to current session folder
    
    % Manual line bisection
    mlbName = dir([dirSess 'MLB_*.mat']); %getting file details for MLB data
    load(mlbName.name);    
    mlb.(sprintf('%s', session)) = data;
    
    % Landmarks
    lmName = dir([dirSess 'LM_*.mat']); %getting file details for MLB data
    load(lmName.name);    
    lm.(sprintf('%s', session)) = data;
    
    % Landmarks 2AFC
    lm2Name= dir([dirSess 'LM2afc_*.mat']); %getting file details for MLB data
    load(lm2Name.name);    
    lm2.(sprintf('%s', session)) = data;
end

%% Analyse MLB data
% Grouping into line length
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));
    % 10 cm line
    mlb.(sprintf('%s', session)).line1mat = ...
        mlb.(sprintf('%s', session)).matrix(find(mlb.(sprintf('%s', session)).matrix(:,1)== 1),:);
    % 20 cm line
    mlb.(sprintf('%s', session)).line2mat = ...
        mlb.(sprintf('%s', session)).matrix(find(mlb.(sprintf('%s', session)).matrix(:,1)== 2),:);
    % 30cm line
    mlb.(sprintf('%s', session)).line3mat = ...
        mlb.(sprintf('%s', session)).matrix(find(mlb.(sprintf('%s', session)).matrix(:,1)== 3),:);
    
    % Average and std error for each line length
    avL1(i) = nanmean(mlb.(sprintf('%s', session)).line1mat(:,6)); %10 cm line
    stdL1(i) = nanstd(mlb.(sprintf('%s', session)).line1mat(:,6));
    mlb.(sprintf('%s', session)).error.line1 = [avL1(i), stdL1(i)];
    avL2(i) = nanmean(mlb.(sprintf('%s', session)).line2mat(:,6)); %20 cm line
    stdL2(i) = nanstd(mlb.(sprintf('%s', session)).line2mat(:,6));
    mlb.(sprintf('%s', session)).error.line2 = [avL2(i), stdL2(i)];
    avL3(i) = nanmean(mlb.(sprintf('%s', session)).line3mat(:,6)); %30 cm line
    stdL3(i) = nanstd(mlb.(sprintf('%s', session)).line3mat(:,6));
    mlb.(sprintf('%s', session)).error.line3 = [avL3(i), stdL3(i)];
end

% Average error across sessions
% For each line length
mlb.line1err = [mean(avL1), std(avL1)];
mlb.line2err = [mean(avL2), std(avL2)];
mlb.line3err = [mean(avL3), std(avL3)];
% Total error and std across sessions
allError = [mlb.line1err(1), mlb.line2err(1), mlb.line3err(1)];
mlb.meanTotError = mean(allError);

%% Analyse LM data 
% Take percentage left-side longer responses for each shift in mm
% Grouping into line length
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));    
    %% Line matrix for each session
    % 10 cm line
    lm1mat = lm.(sprintf('%s', session)).matrix(find(lm.(sprintf('%s', session)).matrix(:,1)== 1),:);
    lm.(sprintf('%s', session)).line1.mat = lm1mat;
    % 20 cm line
    lm2mat = lm.(sprintf('%s', session)).matrix(find(lm.(sprintf('%s', session)).matrix(:,1)== 2),:);
    lm.(sprintf('%s', session)).line2.mat = lm2mat;
    % 30cm line
    lm3mat = lm.(sprintf('%s', session)).matrix(find(lm.(sprintf('%s', session)).matrix(:,1)== 3),:);
    lm.(sprintf('%s', session)).line3.mat = lm3mat;
    
    % Matrix for each shift per line
    lm.(sprintf('%s', session)).line1.mid = lm1mat(find(lm1mat(:,2)== 0),:); %line 1 midpoint
    lm.(sprintf('%s', session)).line2.mid = lm2mat(find(lm2mat(:,2)== 0),:); %line 2 midpoint 
    lm.(sprintf('%s', session)).line3.mid = lm3mat(find(lm3mat(:,2)== 0),:); %line 3 midpoint 
    % For all lines
    sessMat = lm.(sprintf('%s', session)).matrix; %entire session matrix
    lm.(sprintf('%s', session)).allshift.mid = sessMat(find(sessMat(:,2)== 0),:); %midpoint
    for ii = 1:5 %5 shifts per line
        shift = ii*2; shiftmm = shift/10; %for naming
        name = sprintf('%d', shift);
        % Line 1
        lm.(sprintf('%s', session)).line1.(sprintf('left%d', shift))...
            = lm1mat(find(lm1mat(:,2)== 1 & lm1mat(:,3)== shiftmm),:); %left
        lm.(sprintf('%s', session)).line1.(sprintf('right%d', shift))...
            = lm1mat(find(lm1mat(:,2)== 2 & lm1mat(:,3)== shiftmm),:); %right
        % Line 2
        lm.(sprintf('%s', session)).line2.(sprintf('left%d', shift))...
            = lm2mat(find(lm2mat(:,2)== 1 & lm2mat(:,3)== shiftmm),:); %left
        lm.(sprintf('%s', session)).line2.(sprintf('right%d', shift))...
            = lm2mat(find(lm2mat(:,2)== 2 & lm2mat(:,3)== shiftmm),:); %right
        % Line 3
        lm.(sprintf('%s', session)).line3.(sprintf('left%d', shift))...
            = lm3mat(find(lm3mat(:,2)== 1 & lm3mat(:,3)== shiftmm),:); %left
        lm.(sprintf('%s', session)).line3.(sprintf('right%d', shift))...
            = lm3mat(find(lm3mat(:,2)== 2 & lm3mat(:,3)== shiftmm),:); %right
        
        % All lines!
        lm.(sprintf('%s', session)).allshift.(sprintf('left%d', shift))...
            = sessMat(find(sessMat(:,2)== 1 & sessMat(:,3)== shiftmm),:); %left
        lm.(sprintf('%s', session)).allshift.(sprintf('right%d', shift))...
            = sessMat(find(sessMat(:,2)== 2 & sessMat(:,3)== shiftmm),:); %left
    end
    
    %% Find percentage 'left side longer' for each line
    % Vector
    measurements = -10:2:10;
    lm.(sprintf('%s', session)).line1.per(:,1) = measurements'; %line 1
    lm.(sprintf('%s', session)).line2.per(:,1) = measurements'; %line 2
    lm.(sprintf('%s', session)).line3.per(:,1) = measurements'; %line 3
    % Finding the number of responses where left side is longer, summing
    % and then calculating percentage
    perMid = (sum(lm.(sprintf('%s', session)).line1.mid(:,4)==1)/...
        length(lm.(sprintf('%s', session)).line1.mid(:,4)))*100;
    lm.(sprintf('%s', session)).line1.per(1,2) = perMid; %line 1
    perMid = (sum(lm.(sprintf('%s', session)).line2.mid(:,4)==1)/...
        length(lm.(sprintf('%s', session)).line2.mid(:,4)))*100;
    lm.(sprintf('%s', session)).line2.per(1,2) = perMid; %line 2
    perMid = (sum(lm.(sprintf('%s', session)).line3.mid(:,4)==1)/...
        length(lm.(sprintf('%s', session)).line3.mid(:,4)))*100;
    lm.(sprintf('%s', session)).line3.per(1,2) = perMid; %line 3

    % Percentage 'left side longer' for each session (collapsed across
    % line)
    
end



%% Analyse LM data
% Take percentage top-line longer responses for each shift in mm (of the
% top line)
% Grouping into line length
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));
    % 10 cm line
    lm2.(sprintf('%s', session)).line1mat = ...
        lm2.(sprintf('%s', session)).matrix(find(lm2.(sprintf('%s', session)).matrix(:,1)== 1),:);
    % 20 cm line
    lm2.(sprintf('%s', session)).line2mat = ...
        lm2.(sprintf('%s', session)).matrix(find(lm2.(sprintf('%s', session)).matrix(:,1)== 2),:);
    % 30cm line
    lm2.(sprintf('%s', session)).line3mat = ...
        lm2.(sprintf('%s', session)).matrix(find(lm2.(sprintf('%s', session)).matrix(:,1)== 3),:);
end

%% Close and save