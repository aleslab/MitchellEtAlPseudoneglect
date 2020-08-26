%% Analysis for Tactile tasks - AG.Mitchell 01.10.18
% Loads in the data for the task and averages the response across all
% variables
% Then calculates bias
% Then calculates bias relative to line length and shift into L/R hemifield
clear all

%% File paths

%This looks for where the "analyzeData" script is and expects data to be in
%the same directory as that script.  If not it can't find that, it looks
%for the data relative to wherever this script is being run. 
dataLocation = fileparts(which('analyzeData'));
dirData = fullfile(dataLocation,'Data');

if ~exist(dirData,'file')
    disp('cannot find data, trying another path')
    dirData = fullfile(fileparts(mfilename('fullpath')),'..','Data');
    if ~exist(dirData,'file')
        error('Cannot find data')
    end
end


nParticipants = [1:19,21:24,26:30];
%nParticipants = 1;
for p = 1:length(nParticipants)

    ppID = sprintf('P%0*d',2,nParticipants(p)); %for use when navigating files
    matfilename = sprintf('%s_tactileanalysis.mat', ppID);
    nSessions = 1:4; %vector number of sessions each participant does
    
    dirPP = [dirData filesep ppID]; %participant directory
    % Making new anaysis folder for saving
    cd(dirPP)
    mkdir Analysis;
    dirAna = [dirPP filesep 'Analysis' filesep];
    cd(dirAna)
    mkdir Tactile
    dirTact = [dirAna 'Tactile' filesep];

    trb = struct;
    tr2 = struct;

    %% Importing data for both tasks
    for i = 1:length(nSessions)
        session = sprintf('Session%0*d',2,nSessions(i));
        dirSess = [dirPP filesep session filesep]; %current session pathway
        cd(dirSess); %directing to current session folder

        %% TRB task
        [Data,Text] = xlsread('TrialMatrix_TRB'); %importing trial information for landmark tasks
        size = Data(:, ismember(Text(1,:), 'size'));
        side = Data(:, ismember(Text(1,:), 'shift')); %2deg left = 1, middle = 0, 2deg right = 2
        response = Data(:, ismember(Text(1,:), 'response')); %in mm, relative to length of line
        [~, hand] = xlsread('TrialMatrix_TRB','D2:D55'); %hand used L = left, R = right
        % converting hand used from string to numeric, 1= l, 2 = r
        hand_used = [];
        for ii = 1:length(hand)
            if hand{ii} == 'L'
                hand_used(ii,:) = 1;
            elseif hand{ii} == 'R'
                hand_used(ii,:) = 2;
            end
        end
        
        % Finding error 
        size_mm = size*100; %finding actual line length in mm to calculate error
        error = response - size_mm/2; %calculating data error in mm
        error = error/10; %error in cm
        % removing outlier trials from the experiment
        averror = nanmean(error);
        sderror = nanstd(error);
        error(find(error > (averror+2.5*sderror))) = NaN;
        error(find(error < (averror-2.5*sderror))) = NaN;

        trb.(sprintf('%s', session)).data.error = error;
        % data matrix for later analyses
        trb.(sprintf('%s', session)).matrix = [size, side, response, hand_used, error];

        %% TR2 task
        [Data,Text] = xlsread('TrialMatrix_TR2'); %importing trial information for landmark tasks
        size = Data(:, ismember(Text(1,:), 'size'));
        shift = Data(:, ismember(Text(1,:), 'shift')); %shift in mm (0,2,4,6,8,1(=10))
        line = Data(:, ismember(Text(1,:), 'line')); %indicating the line that should be shifted left (the other right, 1 = top, 2 = bottom)
        response = Data(:, ismember(Text(1,:), 'response')); %top line = 1, bottom line = 2
        % Converting shift = 1 to shift = 10, representing actual shift
        shift(find(shift==1)) = 10;
        % Making matrix
        tr2.(sprintf('%s', session)).matrix = [size, shift, line, response];

    end

    %% TRB analysis
    % Grouping into line length
    for i = 1:length(nSessions)
        session = sprintf('Session%0*d',2,nSessions(i));
        cd(dirSess); %directing to current session folder
        length1 = 100; length2 = 200; length3 = 300; %length of lines in mm
        % 10 cm line
        trb.(sprintf('%s', session)).line1mat = ...
            trb.(sprintf('%s', session)).matrix(find(trb.(sprintf('%s', session)).matrix(:,1)== 1),:);
        % 20 cm line
        trb.(sprintf('%s', session)).line2mat = ...
            trb.(sprintf('%s', session)).matrix(find(trb.(sprintf('%s', session)).matrix(:,1)== 2),:);
        % 30cm line
        trb.(sprintf('%s', session)).line3mat = ...
            trb.(sprintf('%s', session)).matrix(find(trb.(sprintf('%s', session)).matrix(:,1)== 3),:);

        % Average and std error for each line length
        avL1(i) = nanmean(trb.(sprintf('%s', session)).line1mat(:,5)); %10 cm line
        stdL1(i) = nanstd(trb.(sprintf('%s', session)).line1mat(:,5));
        trb.(sprintf('%s', session)).error.line1 = [avL1(i), stdL1(i)];
        avL2(i) = nanmean(trb.(sprintf('%s', session)).line2mat(:,5)); %20 cm line
        stdL2(i) = nanstd(trb.(sprintf('%s', session)).line2mat(:,5));
        trb.(sprintf('%s', session)).error.line2 = [avL2(i), stdL2(i)];
        avL3(i) = nanmean(trb.(sprintf('%s', session)).line3mat(:,5)); %30 cm line
        stdL3(i) = nanstd(trb.(sprintf('%s', session)).line3mat(:,5));
        trb.(sprintf('%s', session)).error.line3 = [avL3(i), stdL3(i)];
        
        % Calculating data as proportion of the line (in mm) to standardise
        % for analysis
        propL1 = ((trb.(sprintf('%s', session)).line1mat(:,5))*10)/(length1); %error in mm, as a proportion of half the line length (the error is signed)
        propL2 = ((trb.(sprintf('%s', session)).line2mat(:,5))*10)/(length2);
        propL3 = ((trb.(sprintf('%s', session)).line3mat(:,5))*10)/(length3);
        trb.(sprintf('%s', session)).line1mat(:,6) = propL1; %adding to matrix
        trb.(sprintf('%s', session)).line2mat(:,6) = propL2;
        trb.(sprintf('%s', session)).line3mat(:,6) = propL3;
        
        % Doing the same for hand used
        trb.(sprintf('%s', session)).left_hand = ...
            trb.(sprintf('%s', session)).matrix(find(trb.(sprintf('%s', session)).matrix(:,4)== 1),:); %left hand
        trb.(sprintf('%s', session)).right_hand = ...
            trb.(sprintf('%s', session)).matrix(find(trb.(sprintf('%s', session)).matrix(:,4)== 2),:); %right hand
        % Average and std error for each hand used
        av_left(i) = nanmean(trb.(sprintf('%s', session)).left_hand(:,5)); %10 cm line
        std_left(i) = nanstd(trb.(sprintf('%s', session)).left_hand(:,5));
        trb.(sprintf('%s', session)).error.left_hand = [av_left(i), std_left(i)];
        av_right(i) = nanmean(trb.(sprintf('%s', session)).right_hand(:,5)); %20 cm line
        std_right(i) = nanstd(trb.(sprintf('%s', session)).right_hand(:,5));
        trb.(sprintf('%s', session)).error.right_hand = [av_right(i), std_right(i)];
        
    end

    %% Average error across sessions/lines
    % For each line length
    trb.line1.err = [mean(avL1), std(avL1)];
    trb.line2.err = [mean(avL2), std(avL2)];
    trb.line3.err = [mean(avL3), std(avL3)];
    % Total error and std across sessions
    allError = [trb.line1.err(1), trb.line2.err(1), trb.line3.err(1)];
    trb.meanTotError = mean(allError);
    
    % Matrices for session plotting, each line
    trb.sessionVals.line1.err(1,:) = avL1; trb.sessionVals.line1.err(2,:) = stdL1; %values for sessions
    trb.sessionVals.line2.err(1,:) = avL2; trb.sessionVals.line2.err(2,:) = stdL2;
    trb.sessionVals.line3.err(1,:) = avL3; trb.sessionVals.line3.err(2,:) = stdL3;

    % All sessions
    for i = 1:length(nSessions)
        sess(i,:) = [avL1(1,i), avL2(1,i), avL3(1,i)];
    end
    trb.sessionVals.all.err(1,:) = [mean(sess(1,:)), mean(sess(2,:)),...
        mean(sess(3,:)), mean(sess(4,:))];
    trb.sessionVals.all.err(2,:) = [std(sess(1,:)), std(sess(2,:)),...
        std(sess(3,:)), std(sess(4,:))];
    
    % Same for hand useds
    trb.left_hand.err = [mean(av_left), std(av_left)];
    trb.right_hand.err = [mean(av_right), std(av_right)];
    % Matrices for session plotting, each line
    trb.sessionVals.left_hand.err(1,:) = av_left; trb.sessionVals.left_hand.err(2,:) = std_left; %values for sessions
    trb.sessionVals.right_hand.err(1,:) = av_right; trb.sessionVals.right_hand.err(2,:) = std_right;
    

    %% TRB plots
    cd(dirTact); %navigating to analysis directory to save plots
    for j = 1:3
        line = sprintf('line%d', j);
        % Plotting each line
        figure(j)
        bar(trb.sessionVals.(sprintf('%s', line)).err(1,:))
        hold on
        errorbar(trb.sessionVals.(sprintf('%s', line)).err(1,:), trb.sessionVals.(sprintf('%s', line)).err(2,:), 'k')
        xlabel('Sessions'); ylabel('Bisection error');
        title(sprintf('TRB %s', line));
        saveas(figure(j), sprintf('%s_TRB%s.jpg', ppID, line));
    end
   

    % Plotting average of all lines
    figure(4)
    bar(trb.sessionVals.all.err(1,:))
    hold on
    errorbar(trb.sessionVals.all.err(1,:), trb.sessionVals.all.err(2,:), 'k')
    xlabel('Sessions'); ylabel('Bisection error');
    title('trb all sessions');
    saveas(figure(4), sprintf('%s_TRBsess.jpg', ppID));

    % Plotting averages of all sessions
    lines(1,:) = [trb.line1.err(1), trb.line2.err(1), trb.line3.err(1)];
    lines(2,:) = [trb.line1.err(2), trb.line2.err(2), trb.line3.err(2)];
    figure(5)
    bar(lines(1,:))
    hold on
    errorbar(lines(1,:), lines(2,:), 'k');
    xlabel('Lines'); ylabel('Bisection error');
    title('trb all lines');
    saveas(figure(5), sprintf('%s_TRBlines.jpg', ppID));

    %% TR2 analysis
    %% TR2 plots
    %% Close and save
    close all
    save(matfilename, 'trb', 'tr2');
end

