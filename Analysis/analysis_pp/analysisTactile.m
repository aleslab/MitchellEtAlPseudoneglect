%% Analysis for Tactile tasks - AG.Mitchell 01.10.18
% Loads in the data for the task and averages the response across all
% variables
% Then calculates bias
% Then calculates bias relative to line length and shift into L/R hemifield
clear all

%% File paths
nParticipants = [1:19,21,22,24];
for p = 1:length(nParticipants)
    ppID = sprintf('P%0*d',2,nParticipants(p)); %for use when navigating files
    matfilename = sprintf('%s_tactileanalysis.mat', ppID);
    nSessions = 1:4; %vector number of sessions each participant does
    % Directory
    dirBias = ('M:\Experiments\Bias'); %subject to change depending on where you analyse
    dirPP = [dirBias filesep ppID]; %participant directory
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
        trb.(sprintf('%s', session)).matrix = [size, side, response, error];

        %% TR2 task
        [Data,Text] = xlsread('TrialMatrix_TR2'); %importing trial information for landmark tasks
        size = Data(:, ismember(Text(1,:), 'size'));
        shift = Data(:, ismember(Text(1,:), 'shift')); %shift in mm (0,2,4,6,8,1(=10))
        line = Data(:, ismember(Text(1,:), 'line')); %indicating the line that should be shifted left (the other right, 1 = top, 2 = bottom)
        response = Data(:, ismember(Text(1,:), 'response')); %top line = 1, bottom line = 2
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
        avL1(i) = nanmean(trb.(sprintf('%s', session)).line1mat(:,4)); %10 cm line
        stdL1(i) = nanstd(trb.(sprintf('%s', session)).line1mat(:,4));
        trb.(sprintf('%s', session)).error.line1 = [avL1(i), stdL1(i)];
        avL2(i) = nanmean(trb.(sprintf('%s', session)).line2mat(:,4)); %20 cm line
        stdL2(i) = nanstd(trb.(sprintf('%s', session)).line2mat(:,4));
        trb.(sprintf('%s', session)).error.line2 = [avL2(i), stdL2(i)];
        avL3(i) = nanmean(trb.(sprintf('%s', session)).line3mat(:,4)); %30 cm line
        stdL3(i) = nanstd(trb.(sprintf('%s', session)).line3mat(:,4));
        trb.(sprintf('%s', session)).error.line3 = [avL3(i), stdL3(i)];
        
        % Calculating data as proportion of the line (in mm) to standardise
        % for analysis
        propL1 = ((trb.(sprintf('%s', session)).line1mat(:,4))*10)/(length1); %error in mm, as a proportion of half the line length (the error is signed)
        propL2 = ((trb.(sprintf('%s', session)).line2mat(:,4))*10)/(length2);
        propL3 = ((trb.(sprintf('%s', session)).line3mat(:,4))*10)/(length3);
        trb.(sprintf('%s', session)).line1mat(:,5) = propL1; %adding to matrix
        trb.(sprintf('%s', session)).line2mat(:,5) = propL2;
        trb.(sprintf('%s', session)).line3mat(:,5) = propL3;
        
        % Taking the mean and std of the proportion error
        avpropL1(i) = nanmean(trb.(sprintf('%s', session)).line1mat(:,5)); %10 cm line
        stdpropL1(i) = nanstd(trb.(sprintf('%s', session)).line1mat(:,5));
        trb.(sprintf('%s', session)).proportionError.line1 = [avpropL1(i), stdpropL1(i)];
        avpropL2(i) = nanmean(trb.(sprintf('%s', session)).line2mat(:,5)); %20 cm line
        stdpropL2(i) = nanstd(trb.(sprintf('%s', session)).line2mat(:,5));
        trb.(sprintf('%s', session)).proportionError.line2 = [avpropL2(i), stdpropL2(i)];
        avpropL3(i) = nanmean(trb.(sprintf('%s', session)).line3mat(:,5)); %30 cm line
        stdpropL3(i) = nanstd(trb.(sprintf('%s', session)).line3mat(:,5));
        trb.(sprintf('%s', session)).proportionError.line3 = [avpropL3(i), stdpropL3(i)];     
    end

    %% Average error across sessions/lines
    % For each line length
    trb.line1.err = [mean(avL1), std(avL1)];
    trb.line2.err = [mean(avL2), std(avL2)];
    trb.line3.err = [mean(avL3), std(avL3)];
    % Total error and std across sessions
    allError = [trb.line1.err(1), trb.line2.err(1), trb.line3.err(1)];
    trb.meanTotError = mean(allError);
    % Proportion of line error
    trb.line1.properr = [mean(avpropL1), std(avpropL1)];
    trb.line2.properr = [mean(avpropL2), std(avpropL2)];
    trb.line3.properr = [mean(avpropL3), std(avpropL3)];
    % Total proportion error and std across sessions
    allPropErr = [trb.line1.properr(1), trb.line2.properr(1), trb.line3.properr(1)];
    trb.meanTotPropError = mean(allError);
    
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
    
    % Doing the same for proportion error across all sessions
    trb.sessionVals.line1.properr(1,:) = avpropL1; trb.sessionVals.line1.properr(2,:) = stdpropL1; %values for sessions
    trb.sessionVals.line2.properr(1,:) = avpropL2; trb.sessionVals.line2.properr(2,:) = stdpropL2;
    trb.sessionVals.line3.properr(1,:) = avpropL3; trb.sessionVals.line3.properr(2,:) = stdpropL3;
    % All sessions
    for i = 1:length(nSessions)
        sess(i,:) = [avpropL1(1,i), avpropL2(1,i), avpropL3(1,i)];
    end
    trb.sessionVals.all.properr(1,:) = [mean(sess(1,:)), mean(sess(2,:)),...
        mean(sess(3,:)), mean(sess(4,:))];
    trb.sessionVals.all.properr(2,:) = [std(sess(1,:)), std(sess(2,:)),...
        std(sess(3,:)), std(sess(4,:))];

    %% TRB plots
    cd(dirTact); %navigating to analysis directory to save plots
    for j = 1:3
        line = sprintf('line%d', j);
        % Plotting each line
        figure(j)
        bar(trb.sessionVals.(sprintf('%s', line)).properr(1,:))
        hold on
        errorbar(trb.sessionVals.(sprintf('%s', line)).properr(1,:), trb.sessionVals.(sprintf('%s', line)).properr(2,:), 'k')
        xlabel('Sessions'); ylabel('Proportion bisection error');
        title(sprintf('TRB %s', line));
        saveas(figure(j), sprintf('%s_TRB%s.jpg', ppID, line));
    end
   

    % Plotting average of all lines
    figure(4)
    bar(trb.sessionVals.all.properr(1,:))
    hold on
    errorbar(trb.sessionVals.all.properr(1,:), trb.sessionVals.all.properr(2,:), 'k')
    ylim([-0.3 0.3]);
    xlabel('Sessions'); ylabel('Proportion bisection error');
    title('trb all sessions');
    saveas(figure(4), sprintf('%s_TRBsess.jpg', ppID));

    % Plotting averages of all sessions
    lines(1,:) = [trb.line1.properr(1), trb.line2.properr(1), trb.line3.properr(1)];
    lines(2,:) = [trb.line1.properr(2), trb.line2.properr(2), trb.line3.properr(2)];
    figure(5)
    bar(lines(1,:))
    hold on
    errorbar(lines(1,:), lines(2,:), 'k');
    ylim([-0.3 0.3]);
    xlabel('Lines'); ylabel('Proportion bisection error');
    title('trb all lines');
    saveas(figure(5), sprintf('%s_TRBlines.jpg', ppID));

    %% TR2 analysis
    %% TR2 plots
    %% Close and save
    close all
    save(matfilename, 'trb', 'tr2');
end

