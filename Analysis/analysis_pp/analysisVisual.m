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
    lm2.(sprintf('%s', session)).response = response.key; %adding actual responses to aid with analysis
    lm2.(sprintf('%s', session)).catch = stim.side; %adding stimulus side for lapse trial information
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

%% Average error across sessions
% For each line length
mlb.line1err = [mean(avL1), std(avL1)];
mlb.line2err = [mean(avL2), std(avL2)];
mlb.line3err = [mean(avL3), std(avL3)];
% Total error and std across sessions
allError = [mlb.line1err(1), mlb.line2err(1), mlb.line3err(1)];
mlb.meanTotError = mean(allError);

%% Plotting MLB task

%% Analyse LM data 
% Take percentage left-side longer responses for each shift in mm
% Grouping into line length
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));    
    %% Line matrix for each session
    % 10 cm line
    lmmat1 = lm.(sprintf('%s', session)).matrix(find(lm.(sprintf('%s', session)).matrix(:,1)== 1),:);
    lm.(sprintf('%s', session)).line1.mat = lmmat1;
    % 20 cm line
    lmmat2 = lm.(sprintf('%s', session)).matrix(find(lm.(sprintf('%s', session)).matrix(:,1)== 2),:);
    lm.(sprintf('%s', session)).line2.mat = lmmat2;
    % 30cm line
    lmmat3 = lm.(sprintf('%s', session)).matrix(find(lm.(sprintf('%s', session)).matrix(:,1)== 3),:);
    lm.(sprintf('%s', session)).line3.mat = lmmat3;
    
    % Matrix for each shift per line
    lm.(sprintf('%s', session)).line1.mid = lmmat1(find(lmmat1(:,2)== 0),:); %line 1 midpoint
    lm.(sprintf('%s', session)).line2.mid = lmmat2(find(lmmat2(:,2)== 0),:); %line 2 midpoint 
    lm.(sprintf('%s', session)).line3.mid = lmmat3(find(lmmat3(:,2)== 0),:); %line 3 midpoint 
    % For all lines
    sessMat = lm.(sprintf('%s', session)).matrix; %entire session matrix
    lm.(sprintf('%s', session)).allshift.mid = sessMat(find(sessMat(:,2)== 0),:); %midpoint
    for ii = 1:5 %5 shifts per line
        shift = ii*2; shiftmm = shift/10; %for naming
        name = sprintf('%d', shift);
        % Line 1
        lm.(sprintf('%s', session)).line1.(sprintf('left%d', shift))...
            = lmmat1(find(lmmat1(:,2)== 1 & lmmat1(:,3)== shiftmm),:); %left
        lm.(sprintf('%s', session)).line1.(sprintf('right%d', shift))...
            = lmmat1(find(lmmat1(:,2)== 2 & lmmat1(:,3)== shiftmm),:); %right
        % Line 2
        lm.(sprintf('%s', session)).line2.(sprintf('left%d', shift))...
            = lmmat2(find(lmmat2(:,2)== 1 & lmmat2(:,3)== shiftmm),:); %left
        lm.(sprintf('%s', session)).line2.(sprintf('right%d', shift))...
            = lmmat2(find(lmmat2(:,2)== 2 & lmmat2(:,3)== shiftmm),:); %right
        % Line 3
        lm.(sprintf('%s', session)).line3.(sprintf('left%d', shift))...
            = lmmat3(find(lmmat3(:,2)== 1 & lmmat3(:,3)== shiftmm),:); %left
        lm.(sprintf('%s', session)).line3.(sprintf('right%d', shift))...
            = lmmat3(find(lmmat3(:,2)== 2 & lmmat3(:,3)== shiftmm),:); %right
        
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
    lm.(sprintf('%s', session)).line1.per(6,2) = perMid; %line 1
    perMid = (sum(lm.(sprintf('%s', session)).line2.mid(:,4)==1)/...
        length(lm.(sprintf('%s', session)).line2.mid(:,4)))*100;
    lm.(sprintf('%s', session)).line2.per(6,2) = perMid; %line 2
    perMid = (sum(lm.(sprintf('%s', session)).line3.mid(:,4)==1)/...
        length(lm.(sprintf('%s', session)).line3.mid(:,4)))*100;
    lm.(sprintf('%s', session)).line3.per(6,2) = perMid; %line 3
    
    % Percentage for offsets
    % Vector of values to add left shift into correct location
    k = 1:5; k = flipud(k');
    for ii = 1:5
        shift = ii*2; shiftmm = shift/10; %for naming
        name = sprintf('%d', shift);
        % Line 1
        perL = (sum(lm.(sprintf('%s', session)).line1.(sprintf('left%d', shift))(:,4)==1)/...
        length(lm.(sprintf('%s', session)).line1.(sprintf('left%d', shift))(:,4)))*100; %left
        lm.(sprintf('%s', session)).line1.per(k(ii), 2) = perL;
        perR = (sum(lm.(sprintf('%s', session)).line1.(sprintf('right%d', shift))(:,4)==1)/...
        length(lm.(sprintf('%s', session)).line1.(sprintf('right%d', shift))(:,4)))*100; %right
        lm.(sprintf('%s', session)).line1.per(ii+6, 2) = perR;
        
        %Line 2
        perL = (sum(lm.(sprintf('%s', session)).line2.(sprintf('left%d', shift))(:,4)==1)/...
        length(lm.(sprintf('%s', session)).line2.(sprintf('left%d', shift))(:,4)))*100; %left
        lm.(sprintf('%s', session)).line2.per(k(ii), 2) = perL;
        perR = (sum(lm.(sprintf('%s', session)).line2.(sprintf('right%d', shift))(:,4)==1)/...
        length(lm.(sprintf('%s', session)).line2.(sprintf('right%d', shift))(:,4)))*100; %right
        lm.(sprintf('%s', session)).line2.per(ii+6, 2) = perR;
        
        %Line 3
        perL = (sum(lm.(sprintf('%s', session)).line3.(sprintf('left%d', shift))(:,4)==1)/...
        length(lm.(sprintf('%s', session)).line3.(sprintf('left%d', shift))(:,4)))*100; %left
        lm.(sprintf('%s', session)).line3.per(k(ii), 2) = perL;
        perR = (sum(lm.(sprintf('%s', session)).line3.(sprintf('right%d', shift))(:,4)==1)/...
        length(lm.(sprintf('%s', session)).line3.(sprintf('right%d', shift))(:,4)))*100; %right
        lm.(sprintf('%s', session)).line3.per(ii+6, 2) = perR;       
    end

    %% Percentage 'left side longer' for each session (collapsed across line)
    lm.(sprintf('%s', session)).per(:,1) = measurements';
    % Percenttrage for mid bisected trials for all lines
    perSmid = (sum(lm.(sprintf('%s', session)).allshift.mid(:,4)==1)/...
        length(lm.(sprintf('%s', session)).allshift.mid(:,4)))*100;
    lm.(sprintf('%s', session)).per(6,2) = perSmid;
    
    % For all shifts
    for ii = 1:5
        shift = ii*2; shiftmm = shift/10; %for naming
        name = sprintf('%d', shift);
        
        perSL = (sum(lm.(sprintf('%s', session)).allshift.(sprintf('left%d', shift))(:,4)==1)/...
            length(lm.(sprintf('%s', session)).allshift.(sprintf('left%d', shift))(:,4)))*100; %left
        lm.(sprintf('%s', session)).per(k(ii),2) = perSL;
        perSR = (sum(lm.(sprintf('%s', session)).allshift.(sprintf('right%d', shift))(:,4)==1)/...
            length(lm.(sprintf('%s', session)).allshift.(sprintf('right%d', shift))(:,4)))*100; %right
        lm.(sprintf('%s', session)).per(ii+6,2) = perSR;        
    end
end
% Data across all sessions
% Matrices of shifts for each line length
% Middle bisection
lm.allsessions.line1.mid = [lm.Session01.line1.mid; lm.Session02.line1.mid;...
    lm.Session03.line1.mid; lm.Session04.line1.mid]; %line1
lm.allsessions.line2.mid = [lm.Session01.line2.mid; lm.Session02.line2.mid;...
    lm.Session03.line2.mid; lm.Session04.line2.mid]; %line2
lm.allsessions.line3.mid = [lm.Session01.line3.mid; lm.Session02.line3.mid;...
    lm.Session03.line3.mid; lm.Session04.line3.mid]; %line3
lm.allsessions.alllines.mid = [lm.allsessions.line1.mid; lm.allsessions.line2.mid; lm.allsessions.line3.mid]; %all lines

% Percentage for the middle
% Line 1
lm.allsessions.line1.per(:,1) = measurements';
lm.allsessions.line1.per(6,2) = (sum(lm.allsessions.line1.mid(:,4)==1)/...
    length(lm.allsessions.line1.mid(:,4)))*100;
% Line 2
lm.allsessions.line2.per(:,1) = measurements';
lm.allsessions.line2.per(6,2) = (sum(lm.allsessions.line2.mid(:,4)==1)/...
    length(lm.allsessions.line2.mid(:,4)))*100;
% Line 3
lm.allsessions.line3.per(:,1) = measurements';
lm.allsessions.line3.per(6,2) = (sum(lm.allsessions.line3.mid(:,4)==1)/...
    length(lm.allsessions.line3.mid(:,4)))*100;
% All
lm.allsessions.per(:,1) = measurements';
lm.allsessions.per(6,2) = (sum(lm.allsessions.alllines.mid(:,4)==1)/...
    length(lm.allsessions.alllines.mid(:,4)))*100;

% All other offsets
for ii = 1:5 %number of offsets for left/right
    shift = ii*2; shiftmm = shift/10; %for naming
    name = sprintf('%d', shift);
    % Line 1
    % Matrices
    lm.allsessions.line1.(sprintf('left%d', shift)) = [lm.Session01.line1.(sprintf('left%d', shift)); ...
        lm.Session02.line1.(sprintf('left%d', shift)); lm.Session03.line1.(sprintf('left%d', shift)); ...
        lm.Session04.line1.(sprintf('left%d', shift))]; %left
    lm.allsessions.line1.(sprintf('right%d', shift)) = [lm.Session01.line1.(sprintf('right%d', shift)); ...
        lm.Session02.line1.(sprintf('right%d', shift)); lm.Session03.line1.(sprintf('right%d', shift)); ...
        lm.Session04.line1.(sprintf('right%d', shift))]; %right
    % Percentages
    lm.allsessions.line1.per(k(ii),2) = (sum(lm.allsessions.line1.(sprintf('left%d', shift))(:,4)==1)/...
        length(lm.allsessions.line1.(sprintf('left%d', shift))(:,4)))*100; %left
    lm.allsessions.line1.per(ii+6,2) = (sum(lm.allsessions.line1.(sprintf('right%d', shift))(:,4)==1)/...
        length(lm.allsessions.line1.(sprintf('right%d', shift))(:,4)))*100; %right
    
    % Line 2
    % Matrices
    lm.allsessions.line2.(sprintf('left%d', shift)) = [lm.Session01.line2.(sprintf('left%d', shift)); ...
        lm.Session02.line2.(sprintf('left%d', shift)); lm.Session03.line2.(sprintf('left%d', shift)); ...
        lm.Session04.line2.(sprintf('left%d', shift))]; %left
    lm.allsessions.line2.(sprintf('right%d', shift)) = [lm.Session01.line2.(sprintf('right%d', shift)); ...
        lm.Session02.line2.(sprintf('right%d', shift)); lm.Session03.line2.(sprintf('right%d', shift)); ...
        lm.Session04.line2.(sprintf('right%d', shift))]; %right
    % Percentages
    lm.allsessions.line2.per(k(ii),2) = (sum(lm.allsessions.line2.(sprintf('left%d', shift))(:,4)==1)/...
        length(lm.allsessions.line2.(sprintf('left%d', shift))(:,4)))*100; %left
    lm.allsessions.line2.per(ii+6,2) = (sum(lm.allsessions.line2.(sprintf('right%d', shift))(:,4)==1)/...
        length(lm.allsessions.line2.(sprintf('right%d', shift))(:,4)))*100; %right
    
    % Line 3
    % Matrices
    lm.allsessions.line3.(sprintf('left%d', shift)) = [lm.Session01.line3.(sprintf('left%d', shift)); ...
        lm.Session02.line3.(sprintf('left%d', shift)); lm.Session03.line3.(sprintf('left%d', shift)); ...
        lm.Session04.line3.(sprintf('left%d', shift))]; %left
    lm.allsessions.line3.(sprintf('right%d', shift)) = [lm.Session01.line3.(sprintf('right%d', shift)); ...
        lm.Session02.line3.(sprintf('right%d', shift)); lm.Session03.line3.(sprintf('right%d', shift)); ...
        lm.Session04.line3.(sprintf('right%d', shift))]; %right
    % Percentages
    lm.allsessions.line3.per(k(ii),2) = (sum(lm.allsessions.line3.(sprintf('left%d', shift))(:,4)==1)/...
        length(lm.allsessions.line3.(sprintf('left%d', shift))(:,4)))*100; %left
    lm.allsessions.line3.per(ii+6,2) = (sum(lm.allsessions.line3.(sprintf('right%d', shift))(:,4)==1)/...
        length(lm.allsessions.line3.(sprintf('right%d', shift))(:,4)))*100; %right
    
    % All lines
    % Matrices
    lm.allsessions.alllines.(sprintf('left%d', shift)) = [lm.allsessions.line1.(sprintf('left%d', shift));...
        lm.allsessions.line2.(sprintf('left%d', shift)); lm.allsessions.line3.(sprintf('left%d', shift))]; %left
    lm.allsessions.alllines.(sprintf('right%d', shift)) = [lm.allsessions.line1.(sprintf('right%d', shift));...
        lm.allsessions.line2.(sprintf('right%d', shift)); lm.allsessions.line3.(sprintf('right%d', shift))]; %right
    % Percentages
    lm.allsessions.per(k(ii),2) = (sum(lm.allsessions.alllines.(sprintf('left%d', shift))(:,4)==1)/...
        length(lm.allsessions.alllines.(sprintf('left%d', shift))(:,4)))*100; %left
    lm.allsessions.per(ii+6,2) = (sum(lm.allsessions.alllines.(sprintf('right%d', shift))(:,4)==1)/...
        length(lm.allsessions.alllines.(sprintf('right%d', shift))(:,4)))*100; %left
end

%% Lapse rate for LM task 
% calculating response error for the lines that are clearly longer on l/r
% (10mm shift)
errorL = (sum(lm.allsessions.alllines.left10(:,4)==2)/...
    length(lm.allsessions.alllines.left10(:,4)))*100; %left 10mm shift, finding incorrect responses
errorR = (sum(lm.allsessions.alllines.right10(:,4)==1)/...
    length(lm.allsessions.alllines.right10(:,4)))*100; %right 10mm shift, finding incorrect responses
lm.lapse = mean(errorL, errorR); %lapse rate calculated in percentage

%% Plotting LM task

%% Analyse LM data
% Take percentage top-line longer responses for each shift in mm (of the
% top line)
% Grouping into line length
for i = 1:length(nSessions)
    session = sprintf('Session%0*d',2,nSessions(i));
    % Adding column 8 of matrix to reflect actual response (11 - 1, top; 12
    % - 2, bottom)
    lm2response = lm2.(sprintf('%s', session)).response';
    lm2response(lm2response==11) = 1; lm2response(lm2response==12) = 2; %converting to new values
    lm2.(sprintf('%s', session)).matrix(:,8) = lm2response; %adding to matrix
        
    % First thing's first, need to create separate matrices with lapse rate
    % trials, and then remove them from results matrix
    % Add column to matrix
    lm2.(sprintf('%s', session)).matrix(:,7) = lm2.(sprintf('%s', session)).catch;
    mat = lm2.(sprintf('%s', session)).matrix;
     % Create lapse trial matrix
    lm2.lapse.(sprintf('%sMat', session)) = mat(find(mat(:,7)==1 | mat(:,7)==2),:);
    
    % Getting matrices, and removing lapse trials
    % 10 cm line
    lm2mat1 = lm2.(sprintf('%s', session)).matrix(find(lm2.(sprintf('%s', session)).matrix(:,1)== 1 ...
        & lm2.(sprintf('%s', session)).matrix(:,7)==0),:);
    lm2.(sprintf('%s', session)).line1.mat = lm2mat1;
    % 20 cm line
    lm2mat2 = lm2.(sprintf('%s', session)).matrix(find(lm2.(sprintf('%s', session)).matrix(:,1)== 2 ...
        & lm2.(sprintf('%s', session)).matrix(:,7)==0),:);
    lm2.(sprintf('%s', session)).line2.mat = lm2mat2;
    % 30cm line
    lm2mat3 = lm2.(sprintf('%s', session)).matrix(find(lm2.(sprintf('%s', session)).matrix(:,1)== 3 ...
        & lm2.(sprintf('%s', session)).matrix(:,7)==0),:);
    lm2.(sprintf('%s', session)).line3.mat = lm2mat3;
    
    % Matrix for each shift per line
    % Left-shift line midpoint
    lm2.(sprintf('%s', session)).line1.mid = lm2mat1(find(lm2mat1(:,2)== 0 | ...
        lm2mat1(:,3)== 0),:); %line 1 midpoint
    lm2.(sprintf('%s', session)).line2.mid = lm2mat2(find(lm2mat2(:,2)== 0 |...
        lm2mat2(:,3)== 0),:); %line 2 midpoint 
    lm2.(sprintf('%s', session)).line3.mid = lm2mat3(find(lm2mat3(:,2)== 0 |...
        lm2mat3(:,3)== 0),:); %line 3 midpoint 
    % For all lines
    sessMat = lm2.(sprintf('%s', session)).matrix; %entire session matrix
    sessMat = sessMat(find(sessMat(:,7)==0),:); %removing lapse trials
    lm2.(sprintf('%s', session)).allshift.mid = sessMat(find(sessMat(:,2)== 0 |...
        sessMat(:,3)== 0),:); %midpoint

   for ii = 1:5 %5 shifts per line
        shift = ii*2; shiftmm = shift/10; %for naming
        name = sprintf('%d', shift);
        % Line 1
        lm2.(sprintf('%s', session)).line1.(sprintf('left%d', shift))...
            = lm2mat1(find(lm2mat1(:,2)== shiftmm),:); %top line left
        lm2.(sprintf('%s', session)).line1.(sprintf('right%d', shift))...
            = lm2mat1(find(lm2mat1(:,3)== shiftmm),:); %right
        % Line 2
        lm2.(sprintf('%s', session)).line2.(sprintf('left%d', shift))...
            = lm2mat2(find(lm2mat2(:,2)== shiftmm),:); %left
        lm2.(sprintf('%s', session)).line2.(sprintf('right%d', shift))...
            = lm2mat2(find(lm2mat2(:,3)== shiftmm),:); %right
        % Line 3
        lm2.(sprintf('%s', session)).line3.(sprintf('left%d', shift))...
            = lm2mat3(find(lm2mat3(:,2)== shiftmm),:); %left
        lm2.(sprintf('%s', session)).line3.(sprintf('right%d', shift))...
            = lm2mat3(find(lm2mat3(:,3)== shiftmm),:); %right
        
        % All lines!
        lm2.(sprintf('%s', session)).allshift.(sprintf('left%d', shift))...
            = sessMat(find(sessMat(:,2)== shiftmm),:); %left
        lm2.(sprintf('%s', session)).allshift.(sprintf('right%d', shift))...
            = sessMat(find(sessMat(:,3)== shiftmm),:); %left
   end
    
   %% Find percentage 'left shifted line longer' for each line
    % Vector
    measurements = -10:2:10;
    lm2.(sprintf('%s', session)).line1.per(:,1) = measurements'; %line 1
    lm2.(sprintf('%s', session)).line2.per(:,1) = measurements'; %line 2
    lm2.(sprintf('%s', session)).line3.per(:,1) = measurements'; %line 3
   
    % Calculating vectors of equal values between 'line to shift left' and
    % 'response' - used to be able to calculate percentage of
    % 'left-shitfted line is longer' responses
    % Midpoint lines - for each line
    for j = 1:length(lm2.(sprintf('%s', session)).line1.mid(:,4))
        lm2.(sprintf('%s', session)).line1.mid(j,9) = ...
            isequal(lm2.(sprintf('%s', session)).line1.mid(j,4), lm2.(sprintf('%s', session)).line1.mid(j,8));
        lm2.(sprintf('%s', session)).line2.mid(j,9) = ...
            isequal(lm2.(sprintf('%s', session)).line2.mid(j,4), lm2.(sprintf('%s', session)).line2.mid(j,8));
        lm2.(sprintf('%s', session)).line3.mid(j,9) = ...
            isequal(lm2.(sprintf('%s', session)).line3.mid(j,4), lm2.(sprintf('%s', session)).line3.mid(j,8));
    end
    
    % Calculating the percentage
    perMid = sum((lm2.(sprintf('%s', session)).line1.mid(:,9))/...
        length(lm2.(sprintf('%s', session)).line1.mid(:,9)))*100; %line 1
    lm.(sprintf('%s', session)).line1.per(6, 2) = perMid;
    perMid = sum((lm2.(sprintf('%s', session)).line2.mid(:,9))/...
        length(lm2.(sprintf('%s', session)).line2.mid(:,9)))*100; %line 2
    lm.(sprintf('%s', session)).line2.per(6, 2) = perMid;
    perMid = sum((lm2.(sprintf('%s', session)).line3.mid(:,9))/...
        length(lm2.(sprintf('%s', session)).line3.mid(:,9)))*100; %line 3
    lm.(sprintf('%s', session)).line3.per(6, 2) = perMid;
    
    % Doing the same for the shifted lines
    % Vector of values to add left shift into correct location
    k = 1:5; k = flipud(k');
    for ii = 1:5 %number of shifts
        shift = ii*2; shiftmm = shift/10; %for naming
        name = sprintf('%d', shift);
        
        for j = 1:length(lm2.(sprintf('%s', session)).line1.(sprintf('left%d', shift))(:,4))
            %% continue here...
            % Line 1
            lm2.(sprintf('%s', session)).line1.(sprintf('left%d', shift))(j,9) = ...
                isequal(lm2.(sprintf('%s', session)).line1.(sprintf('left%d', shift))(j,4),...
            lm2.(sprintf('%s', session)).line1.(sprintf('left%d', shift))(j,8)); %left
            lm2.(sprintf('%s', session)).line1.(sprintf('right%d', shift))(j,9) = ...
                isequal(lm2.(sprintf('%s', session)).line1.(sprintf('right%d', shift))(j,4),...
            lm2.(sprintf('%s', session)).line1.(sprintf('right%d', shift))(j,8)); %right
            % Line 2
            lm2.(sprintf('%s', session)).line2.(sprintf('left%d', shift))(j,9) = ...
                isequal(lm2.(sprintf('%s', session)).line2.(sprintf('left%d', shift))(j,4),...
            lm2.(sprintf('%s', session)).line2.(sprintf('left%d', shift))(j,8)); %left
            lm2.(sprintf('%s', session)).line2.(sprintf('right%d', shift))(j,9) = ...
                isequal(lm2.(sprintf('%s', session)).line2.(sprintf('right%d', shift))(j,4),...
            lm2.(sprintf('%s', session)).line2.(sprintf('right%d', shift))(j,8)); %right
            % Line 3
            lm2.(sprintf('%s', session)).line3.(sprintf('left%d', shift))(j,9) = ...
                isequal(lm2.(sprintf('%s', session)).line3.(sprintf('left%d', shift))(j,4),...
            lm2.(sprintf('%s', session)).line3.(sprintf('left%d', shift))(j,8)); %left
            lm2.(sprintf('%s', session)).line3.(sprintf('right%d', shift))(j,9) = ...
                isequal(lm2.(sprintf('%s', session)).line3.(sprintf('right%d', shift))(j,4),...
            lm2.(sprintf('%s', session)).line3.(sprintf('right%d', shift))(j,8)); %right
        end
    end

end
%% Plotting LM2 task
%% Lapse rate for LM2 task
%% Close and save