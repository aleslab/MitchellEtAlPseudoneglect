%% A.G. Mitchell 22.10.18
%% Visual analysis
clear all

%% File paths
nParticipants = [1:19,21:24,26:30];
for p = 1:length(nParticipants)
    ppID = sprintf('P%0*d',2,nParticipants(p)); %for use when navigating files
    %ppID = input('Participant ID? ', 's'); %for use when navigating files, debugging version

    matfilename = sprintf('%s_visualanalysisStart.mat', ppID);
    nSessions = 1:4; %vector number of sessions each participant does
    % Directory
    dirBias = ('M:\Alex_Files\Experiments\Bias'); %subject to change depending on where you analyse
    dirPP = [dirBias filesep 'Data' filesep ppID]; %participant directory
    % Making new anaysis folder for saving
    cd(dirPP)
    mkdir Analysis;
    dirAna = [dirPP filesep 'Analysis' filesep];
    cd(dirAna)
    mkdir Visual
    dirVis = [dirAna 'Visual' filesep];

    % Load data for all trials
    for i = 1:length(nSessions)
        session = sprintf('Session%0*d',2,nSessions(i));
        dirSess = [dirPP filesep session filesep]; %current session pathway
        cd(dirSess); %directing to current session folder

        % Manual line bisection
        mlbName = dir([dirSess 'MLB_*.mat']); %getting file details for MLB data
        load(mlbName.name);    
        mlb.(sprintf('%s', session)) = data;
        hand_used = zeros(length(data.matrix),1);
        
        % Adding 'hand used' information to address reviwer comments
        % (30.04.20) - First 15 participants used right hand, then left;
        % final 15 used left then right
        if nParticipants(p) > 0 && nParticipants(p) < 16
            hand_used(1:90,1) = 2; %right hand used first 
            hand_used(91:end,1) = 1; %left hand second
            mlb.(sprintf('%s', session)).matrix(:,8) = hand_used;
        else
            hand_used(1:90,1) = 1; %right hand used first 
            hand_used(91:end,1) = 2; %left hand second
            mlb.(sprintf('%s', session)).matrix(:,8) = hand_used;
        end

        % Landmarks
        lmName = dir([dirSess 'LM_*.mat']); %getting file details for MLB data
        load(lmName.name);    
        lm.(sprintf('%s', session)) = data;
        lm.(sprintf('%s', session)).matrix(:,5) = stim.pos;
        
        % Landmarks 2AFC
        lm2Name= dir([dirSess 'LM2afc_*.mat']); %getting file details for MLB data
        load(lm2Name.name);    
        lm2.(sprintf('%s', session)) = data;
        lm2.(sprintf('%s', session)).response = response.key; %adding actual responses to aid with analysis
        lm2.(sprintf('%s', session)).catch = stim.side; %adding stimulus side for lapse trial information
    end

    %% Analyse MLB data
    
    for i = 1:length(nSessions)
        % Grouping into line length
        session = sprintf('Session%0*d',2,nSessions(i));
        cd(dirSess); %directing to current session folder
        length1 = 100; length2 = 200; length3 = 300; %length of lines in mm
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
        
        % Calculating data as proportion of the line (in mm) to standardise
        % for analysis
        propL1 = ((mlb.(sprintf('%s', session)).line1mat(:,6))*10)/(length1); %error in mm, as a proportion of half the line length (the error is signed)
        propL2 = ((mlb.(sprintf('%s', session)).line2mat(:,6))*10)/(length2);
        propL3 = ((mlb.(sprintf('%s', session)).line3mat(:,6))*10)/(length3);
        mlb.(sprintf('%s', session)).line1mat(:,8) = propL1; %adding to matrix
        mlb.(sprintf('%s', session)).line2mat(:,8) = propL2;
        mlb.(sprintf('%s', session)).line3mat(:,8) = propL3;
        
        % Taking the mean and std of the proportion error
        avpropL1(i) = nanmean(mlb.(sprintf('%s', session)).line1mat(:,8)); %10 cm line
        stdpropL1(i) = nanstd(mlb.(sprintf('%s', session)).line1mat(:,8));
        mlb.(sprintf('%s', session)).proportionError.line1 = [avpropL1(i), stdpropL1(i)];
        avpropL2(i) = nanmean(mlb.(sprintf('%s', session)).line2mat(:,8)); %20 cm line
        stdpropL2(i) = nanstd(mlb.(sprintf('%s', session)).line2mat(:,8));
        mlb.(sprintf('%s', session)).proportionError.line2 = [avpropL2(i), stdpropL2(i)];
        avpropL3(i) = nanmean(mlb.(sprintf('%s', session)).line3mat(:,8)); %30 cm line
        stdpropL3(i) = nanstd(mlb.(sprintf('%s', session)).line3mat(:,8));
        mlb.(sprintf('%s', session)).proportionError.line3 = [avpropL3(i), stdpropL3(i)]; 
        
        % grouping into hand used
        % left hand
        mlb.(sprintf('%s', session)).left_hand = ...
            mlb.(sprintf('%s', session)).matrix(find(mlb.(sprintf('%s', session)).matrix(:,8)== 1),:);
        % right hand
        mlb.(sprintf('%s', session)).right_hand = ...
            mlb.(sprintf('%s', session)).matrix(find(mlb.(sprintf('%s', session)).matrix(:,8)== 2),:);
        % Average and std error for each hand used
        av_leftHand(i) = nanmean(mlb.(sprintf('%s', session)).left_hand(:,6)); %left hand
        std_leftHand(i) = nanstd(mlb.(sprintf('%s', session)).left_hand(:,6));
        mlb.(sprintf('%s', session)).error.left_hand = [av_leftHand(i), std_leftHand(i)];
        av_rightHand(i) = nanmean(mlb.(sprintf('%s', session)).right_hand(:,6)); %right hand
        std_rightHand(i) = nanstd(mlb.(sprintf('%s', session)).right_hand(:,6));
        mlb.(sprintf('%s', session)).error.right_hand = [av_rightHand(i), std_rightHand(i)];
    end

    %% Average error across sessions
    % For each line length
    mlb.line1.err = [mean(avL1), std(avL1)];
    mlb.line2.err = [mean(avL2), std(avL2)];
    mlb.line3.err = [mean(avL3), std(avL3)];
    % Total error and std across sessions
    allError = [mlb.line1.err(1), mlb.line2.err(1), mlb.line3.err(1)];
    mlb.meanTotError = mean(allError);
    % Proportion of line error
    mlb.line1.properr = [mean(avpropL1), std(avpropL1)];
    mlb.line2.properr = [mean(avpropL2), std(avpropL2)];
    mlb.line3.properr = [mean(avpropL3), std(avpropL3)];
    % Total proportion error and std across sessions
    allPropErr = [mlb.line1.properr(1), mlb.line2.properr(1), mlb.line3.properr(1)];
    mlb.meanTotPropError = mean(allError);
    
    % The same for hand used 
    mlb.left_hand.err = [mean(av_leftHand), std(av_leftHand)];
    mlb.right_hand.err = [mean(av_rightHand), std(av_rightHand)];
    
    %% Plotting MLB task
    cd(dirVis); %navigating to analysis folder to save plots to
    % Plot across all sessions
    % Matrices for session plotting, each line
    mlb.sessionVals.line1.err(1,:) = avL1; mlb.sessionVals.line1.err(2,:) = stdL1; %values for sessions
    mlb.sessionVals.line2.err(1,:) = avL2; mlb.sessionVals.line2.err(2,:) = stdL2;
    mlb.sessionVals.line3.err(1,:) = avL3; mlb.sessionVals.line3.err(2,:) = stdL3;
    % All sessions
    for i = 1:length(nSessions)
        sess(i,:) = [avL1(1,i), avL2(1,i), avL3(1,i)];
    end
    mlb.sessionVals.all.err(1,:) = [mean(sess(1,:)), mean(sess(2,:)),...
        mean(sess(3,:)), mean(sess(4,:))];
    mlb.sessionVals.all.err(2,:) = [std(sess(1,:)), std(sess(2,:)),...
        std(sess(3,:)), std(sess(4,:))];
    
    % Doing the same for proportion error across all sessions
    mlb.sessionVals.line1.properr(1,:) = avpropL1; mlb.sessionVals.line1.properr(2,:) = stdpropL1; %values for sessions
    mlb.sessionVals.line2.properr(1,:) = avpropL2; mlb.sessionVals.line2.properr(2,:) = stdpropL2;
    mlb.sessionVals.line3.properr(1,:) = avpropL3; mlb.sessionVals.line3.properr(2,:) = stdpropL3;
    % All sessions
    for i = 1:length(nSessions)
        sess(i,:) = [avpropL1(1,i), avpropL2(1,i), avpropL3(1,i)];
    end
    mlb.sessionVals.all.properr(1,:) = [mean(sess(1,:)), mean(sess(2,:)),...
        mean(sess(3,:)), mean(sess(4,:))];
    mlb.sessionVals.all.properr(2,:) = [std(sess(1,:)), std(sess(2,:)),...
        std(sess(3,:)), std(sess(4,:))];

    for j = 1:3
        line = sprintf('line%d', j);
        % Plotting each line
        figure(j)
        bar(mlb.sessionVals.(sprintf('%s', line)).properr(1,:))
        hold on
        errorbar(mlb.sessionVals.(sprintf('%s', line)).properr(1,:), mlb.sessionVals.(sprintf('%s', line)).properr(2,:), 'k')
        xlabel('Sessions'); ylabel('Proportion bisection error');
        title(sprintf('MLB %s', line));
        saveas(figure(j), sprintf('%s_MLB%s.jpg', ppID, line));
    end

    % Plotting average of all lines
    figure(4)
    bar(mlb.sessionVals.all.properr(1,:))
    hold on
    errorbar(mlb.sessionVals.all.properr(1,:), mlb.sessionVals.all.properr(2,:), 'k')
    ylim([-0.2 0.2]);
    xlabel('Sessions'); ylabel('Proportion bisection error');
    title('MLB all sessions');
    saveas(figure(4), sprintf('%s_MLBsess.jpg', ppID));

    % Plotting averages of all sessions
    lines(1,:) = [mlb.line1.properr(1), mlb.line2.properr(1), mlb.line3.properr(1)];
    lines(2,:) = [mlb.line1.properr(2), mlb.line2.properr(2), mlb.line3.properr(2)];
    figure(5)
    bar(lines(1,:))
    hold on
    errorbar(lines(1,:), lines(2,:), 'k');
    ylim([-0.2 0.2]);
    xlabel('Lines'); ylabel('Bisection error');
    title('MLB all lines');
    saveas(figure(5), sprintf('%s_MLBlines.jpg', ppID));


    %% Analyse LM data 
    % Take percentage left-side longer responses for each shift in mm
    % Grouping into line length
    for i = 1:length(nSessions)
        session = sprintf('Session%0*d',2,nSessions(i));
        cd(dirSess); %directing to current session folder
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

        %% Find percentage 'right side longer' for each line
        % Structure of results file: [stim aysmmetry, ntrials per cond, n right-side longer per cond, % per cond]
        % Vector
        measurements = -10:2:10; %data to go along the x-axis, stimulus asymmetries
        % Finding the number of responses where left side is longer, summing
        % and then calculating percentage
        k = 1:5; k = flipud(k');
        for j = 1:3 %number of lines
            line = sprintf('line%d', j);
            % getting important variables
            lm.(sprintf('%s', session)).(sprintf('%s', line)).res(:,1) = measurements'; %stimulus asymmetry
            % Midpoint
            nTrials = length(lm.(sprintf('%s', session)).(sprintf('%s', line)).mid(:,4)); %number of trials for  the condition
            nRight = sum(lm.(sprintf('%s', session)).(sprintf('%s', line)).mid(:,4)==2); %number of trials they responded that the right side is longer
            percentR = (nRight/nTrials)*100;
            % adding them to key matrix (just middle)
            lm.(sprintf('%s', session)).(sprintf('%s', line)).res(6,2) = nTrials;
            lm.(sprintf('%s', session)).(sprintf('%s', line)).res(6,3) = nRight;
            lm.(sprintf('%s', session)).(sprintf('%s', line)).res(6,4) = percentR;

            % For L/R shifts
            for ii = 1:5
                shift = ii*2; shiftmm = shift/10; %for naming
                name = sprintf('%d', shift);
                % Getting key values, left shift
                nTrialsL = length(lm.(sprintf('%s', session)).(sprintf('%s', line)).(sprintf('left%d', shift))(:,4)); %number of trials for  the condition
                nRightL = sum(lm.(sprintf('%s', session)).(sprintf('%s', line)).(sprintf('left%d', shift))(:,4)==2); %number of trials they responded that the right side is longer
                percentRL = (nRightL/nTrialsL)*100;
                % Getting key values, right shift
                nTrialsR = length(lm.(sprintf('%s', session)).(sprintf('%s', line)).(sprintf('right%d', shift))(:,4)); %number of trials for  the condition
                nRightR = sum(lm.(sprintf('%s', session)).(sprintf('%s', line)).(sprintf('right%d', shift))(:,4)==2); %number of trials they responded that the right side is longer
                percentRR = (nRightR/nTrialsR)*100;
                % Adding to matrix, left shift
                lm.(sprintf('%s', session)).(sprintf('%s', line)).res(k(ii),2) = nTrialsL;
                lm.(sprintf('%s', session)).(sprintf('%s', line)).res(k(ii),3) = nRightL;
            	lm.(sprintf('%s', session)).(sprintf('%s', line)).res(k(ii),4) = percentRL;
                % Adding to matrix, left shift
                lm.(sprintf('%s', session)).(sprintf('%s', line)).res(ii+6,2) = nTrialsR;
                lm.(sprintf('%s', session)).(sprintf('%s', line)).res(ii+6,3) = nRightR;
            	lm.(sprintf('%s', session)).(sprintf('%s', line)).res(ii+6,4) = percentRR;    
            end
        end

        %% Percentage 'right side longer' for each session (collapsed across line)
        lm.(sprintf('%s', session)).res(:,1) = measurements';
        nTrials = length(lm.(sprintf('%s', session)).allshift.mid(:,4)); %number of trials for  the condition
        nRight = sum(lm.(sprintf('%s', session)).allshift.mid(:,4)==2); %number of trials they responded that the right side is longer
        percentR = (nRight/nTrials)*100;
        % Adding to matrix
        lm.(sprintf('%s', session)).res(6,2) = nTrials;
        lm.(sprintf('%s', session)).res(6,3) = nRight;
        lm.(sprintf('%s', session)).res(6,4) = percentR;

        % For all shifts
        for ii = 1:5
            shift = ii*2; shiftmm = shift/10; %for naming
            name = sprintf('%d', shift);

            % Getting key values, left shift
            nTrialsL = length(lm.(sprintf('%s', session)).allshift.(sprintf('left%d', shift))(:,4)); %number of trials for  the condition
            nRightL = sum(lm.(sprintf('%s', session)).allshift.(sprintf('left%d', shift))(:,4)==2); %number of trials they responded that the right side is longer
            percentRL = (nRightL/nTrialsL)*100;
            % Getting key values, right shift
            nTrialsR = length(lm.(sprintf('%s', session)).allshift.(sprintf('right%d', shift))(:,4)); %number of trials for  the condition
            nRightR = sum(lm.(sprintf('%s', session)).allshift.(sprintf('right%d', shift))(:,4)==2); %number of trials they responded that the right side is longer
            percentRR = (nRightR/nTrialsR)*100;
            % Adding to matrix, left shift
            lm.(sprintf('%s', session)).res(k(ii),2) = nTrialsL;
            lm.(sprintf('%s', session)).res(k(ii),3) = nRightL;
            lm.(sprintf('%s', session)).res(k(ii),4) = percentRL;
            % Adding to matrix, left shift
            lm.(sprintf('%s', session)).res(ii+6,2) = nTrialsR;
            lm.(sprintf('%s', session)).res(ii+6,3) = nRightR;
            lm.(sprintf('%s', session)).res(ii+6,4) = percentRR;      
        end
        %% Calculate percentage 'right side longer' for each colour shift (contrast switch at line midpoint)
        % Response to reviewer comment
        % First, extract shift matrix for each contrast
        % isolate each shift
        mid = lm.(sprintf('%s', session)).allshift.mid;
        left2 = lm.(sprintf('%s', session)).allshift.left2;
        left4 = lm.(sprintf('%s', session)).allshift.left4;
        left6 = lm.(sprintf('%s', session)).allshift.left6;
        left8 = lm.(sprintf('%s', session)).allshift.left8;
        left10 = lm.(sprintf('%s', session)).allshift.left10;
        right2 = lm.(sprintf('%s', session)).allshift.right2;
        right4 = lm.(sprintf('%s', session)).allshift.right4;
        right6 = lm.(sprintf('%s', session)).allshift.right6;
        right8 = lm.(sprintf('%s', session)).allshift.right8;
        right10 = lm.(sprintf('%s', session)).allshift.right10;

        % extract each contrast from each shift
        % black = top left
        lm.(sprintf('%s', session)).con1shift.mid = mid(find(mid(:,5) == 1), :);
        lm.(sprintf('%s', session)).con1shift.left2 = left2(find(left2(:,5) == 1), :);
        lm.(sprintf('%s', session)).con1shift.left4 = left4(find(left4(:,5) == 1), :);
        lm.(sprintf('%s', session)).con1shift.left6 = left6(find(left6(:,5) == 1), :);
        lm.(sprintf('%s', session)).con1shift.left8 = left8(find(left8(:,5) == 1), :);
        lm.(sprintf('%s', session)).con1shift.left10 = left10(find(left10(:,5) == 1), :);
        lm.(sprintf('%s', session)).con1shift.right2 = right2(find(right2(:,5) == 1), :);
        lm.(sprintf('%s', session)).con1shift.right4 = right4(find(right4(:,5) == 1), :);
        lm.(sprintf('%s', session)).con1shift.right6 = right6(find(right6(:,5) == 1), :);
        lm.(sprintf('%s', session)).con1shift.right8 = right8(find(right8(:,5) == 1), :);
        lm.(sprintf('%s', session)).con1shift.right10 = right10(find(right10(:,5) == 1), :);
        % white = top left
        lm.(sprintf('%s', session)).con2shift.mid = mid(find(mid(:,5) == 2), :);
        lm.(sprintf('%s', session)).con2shift.left2 = left2(find(left2(:,5) == 2), :);
        lm.(sprintf('%s', session)).con2shift.left4 = left4(find(left4(:,5) == 2), :);
        lm.(sprintf('%s', session)).con2shift.left6 = left6(find(left6(:,5) == 2), :);
        lm.(sprintf('%s', session)).con2shift.left8 = left8(find(left8(:,5) == 2), :);
        lm.(sprintf('%s', session)).con2shift.left10 = left10(find(left10(:,5) == 2), :);
        lm.(sprintf('%s', session)).con2shift.right2 = right2(find(right2(:,5) == 2), :);
        lm.(sprintf('%s', session)).con2shift.right4 = right4(find(right4(:,5) == 2), :);
        lm.(sprintf('%s', session)).con2shift.right6 = right6(find(right6(:,5) == 2), :);
        lm.(sprintf('%s', session)).con2shift.right8 = right8(find(right8(:,5) == 2), :);
        lm.(sprintf('%s', session)).con2shift.right10 = right10(find(right10(:,5) == 2), :);

        % Get matrices of key values - proportion right for each shift
        % contrast 1 = black top left
        lm.(sprintf('%s', session)).con1(:,1) = measurements'; %black = top left    
        nTrials = length(lm.(sprintf('%s', session)).con1shift.mid(:,4)); %number of trials for  the condition
        nRight = sum(lm.(sprintf('%s', session)).con1shift.mid(:,4)==2); %number of trials they responded that the right side is longer
        percentR = (nRight/nTrials)*100;
        % Adding to matrix
        lm.(sprintf('%s', session)).con1(6,2) = nTrials;
        lm.(sprintf('%s', session)).con1(6,3) = nRight;
        lm.(sprintf('%s', session)).con1(6,4) = percentR;
        % For all shifts
        for ii = 1:5
            shift = ii*2; shiftmm = shift/10; %for naming
            name = sprintf('%d', shift);

            % Getting key values, left shift
            nTrialsL = length(lm.(sprintf('%s', session)).con1shift.(sprintf('left%d', shift))(:,4)); %number of trials for  the condition
            nRightL = sum(lm.(sprintf('%s', session)).con1shift.(sprintf('left%d', shift))(:,4)==2); %number of trials they responded that the right side is longer
            percentRL = (nRightL/nTrialsL)*100;
            % Getting key values, right shift
            nTrialsR = length(lm.(sprintf('%s', session)).con1shift.(sprintf('right%d', shift))(:,4)); %number of trials for  the condition
            nRightR = sum(lm.(sprintf('%s', session)).con1shift.(sprintf('right%d', shift))(:,4)==2); %number of trials they responded that the right side is longer
            percentRR = (nRightR/nTrialsR)*100;
            % Adding to matrix, left shift
            lm.(sprintf('%s', session)).con1(k(ii),2) = nTrialsL;
            lm.(sprintf('%s', session)).con1(k(ii),3) = nRightL;
            lm.(sprintf('%s', session)).con1(k(ii),4) = percentRL;
            % Adding to matrix, left shift
            lm.(sprintf('%s', session)).con1(ii+6,2) = nTrialsR;
            lm.(sprintf('%s', session)).con1(ii+6,3) = nRightR;
            lm.(sprintf('%s', session)).con1(ii+6,4) = percentRR;      
        end
        
        % contrast 2 = white top left
        lm.(sprintf('%s', session)).con2(:,1) = measurements'; %white = top left 
        nTrials = length(lm.(sprintf('%s', session)).con2shift.mid(:,4)); %number of trials for  the condition
        nRight = sum(lm.(sprintf('%s', session)).con2shift.mid(:,4)==2); %number of trials they responded that the right side is longer
        percentR = (nRight/nTrials)*100;
        % Adding to matrix
        lm.(sprintf('%s', session)).con2(6,2) = nTrials;
        lm.(sprintf('%s', session)).con2(6,3) = nRight;
        lm.(sprintf('%s', session)).con2(6,4) = percentR;
        % For all shifts
        for ii = 1:5
            shift = ii*2; shiftmm = shift/10; %for naming
            name = sprintf('%d', shift);

            % Getting key values, left shift
            nTrialsL = length(lm.(sprintf('%s', session)).con2shift.(sprintf('left%d', shift))(:,4)); %number of trials for  the condition
            nRightL = sum(lm.(sprintf('%s', session)).con2shift.(sprintf('left%d', shift))(:,4)==2); %number of trials they responded that the right side is longer
            percentRL = (nRightL/nTrialsL)*100;
            % Getting key values, right shift
            nTrialsR = length(lm.(sprintf('%s', session)).con2shift.(sprintf('right%d', shift))(:,4)); %number of trials for  the condition
            nRightR = sum(lm.(sprintf('%s', session)).con2shift.(sprintf('right%d', shift))(:,4)==2); %number of trials they responded that the right side is longer
            percentRR = (nRightR/nTrialsR)*100;
            % Adding to matrix, left shift
            lm.(sprintf('%s', session)).con2(k(ii),2) = nTrialsL;
            lm.(sprintf('%s', session)).con2(k(ii),3) = nRightL;
            lm.(sprintf('%s', session)).con2(k(ii),4) = percentRL;
            % Adding to matrix, left shift
            lm.(sprintf('%s', session)).con2(ii+6,2) = nTrialsR;
            lm.(sprintf('%s', session)).con2(ii+6,3) = nRightR;
            lm.(sprintf('%s', session)).con2(ii+6,4) = percentRR;      
        end
    end

    
    %% Data across all sessions
    % Matrices of shifts for each line length
    for j = 1:3
        line = sprintf('line%d', j);
        % midpoint
        lm.allsessions.(sprintf('%s', line)).mid = [lm.Session01.(sprintf('%s', line)).mid; ...
            lm.Session02.(sprintf('%s', line)).mid; lm.Session03.(sprintf('%s', line)).mid; ...
            lm.Session04.(sprintf('%s', line)).mid]; % getting all session data for all lines  
        lm.allsessions.(sprintf('%s', line)).res(:,1) = measurements'; 
        % Getting important numbers
        nTrials = length(lm.allsessions.(sprintf('%s', line)).mid(:,4)); %number of trials for  the condition
        nRight = sum(lm.allsessions.(sprintf('%s', line)).mid(:,4)==2); %number of trials they responded that the right side is longer
        percentR = (nRight/nTrials)*100;
        lm.allsessions.(sprintf('%s', line)).res(6,2) = nTrials;
        lm.allsessions.(sprintf('%s', line)).res(6,3) = nRight;
        lm.allsessions.(sprintf('%s', line)).res(6,4) = percentR;
        
        for ii = 1:5
            shift = ii*2; shiftmm = shift/10; %for naming
            name = sprintf('%d', shift);
            % All session line matrices
            lm.allsessions.(sprintf('%s', line)).(sprintf('left%d', shift)) = [lm.Session01.(sprintf('%s', line)).(sprintf('left%d', shift)); ...
                lm.Session02.(sprintf('%s', line)).(sprintf('left%d', shift)); lm.Session03.(sprintf('%s', line)).(sprintf('left%d', shift)); ...
                lm.Session04.line1.(sprintf('left%d', shift))]; %left
            lm.allsessions.(sprintf('%s', line)).(sprintf('right%d', shift)) = [lm.Session01.(sprintf('%s', line)).(sprintf('right%d', shift)); ...
                lm.Session02.(sprintf('%s', line)).(sprintf('right%d', shift)); lm.Session03.(sprintf('%s', line)).(sprintf('right%d', shift)); ...
                lm.Session04.(sprintf('%s', line)).(sprintf('right%d', shift))]; %right
            % Key values, left
            nTrialsL = length(lm.allsessions.(sprintf('%s', line)).(sprintf('left%d', shift))(:,4)); %number of trials for  the condition
            nRightL = sum(lm.allsessions.(sprintf('%s', line)).(sprintf('left%d', shift))(:,4)==2); %number of trials they responded that the right side is longer
            percentRL = (nRightL/nTrialsL)*100;
            % Key values, right
            nTrialsR = length(lm.allsessions.(sprintf('%s', line)).(sprintf('right%d', shift))(:,4)); %number of trials for  the condition
            nRightR = sum(lm.allsessions.(sprintf('%s', line)).(sprintf('right%d', shift))(:,4)==2); %number of trials they responded that the right side is longer
            percentRR = (nRightR/nTrialsR)*100;
            % Matrix left
            lm.allsessions.(sprintf('%s', line)).res(k(ii),2) = nTrialsL;
            lm.allsessions.(sprintf('%s', line)).res(k(ii),3) = nRightL;
            lm.allsessions.(sprintf('%s', line)).res(k(ii),4) = percentRL;
            % Matrix right
            lm.allsessions.(sprintf('%s', line)).res(ii+6,2) = nTrialsR;
            lm.allsessions.(sprintf('%s', line)).res(ii+6,3) = nRightR;
            lm.allsessions.(sprintf('%s', line)).res(ii+6,4) = percentRR;
        end
    end
    
    
%% All sessions, for each contrast
% contrast 1 - black top left
% midpoint 
lm.allsessions.con1shift.mid = [lm.Session01.con1shift.mid; lm.Session02.con1shift.mid; ...
    lm.Session03.con1shift.mid; lm.Session04.con1shift.mid];
lm.allsessions.con1(:,1) = measurements';
nTrials = length(lm.allsessions.con1shift.mid(:,4)); 
nRight = sum(lm.allsessions.con1shift.mid(:,4)==2); %number of trials they responded that the right side is longer
percentR = (nRight/nTrials)*100;
% Adding to matrix
lm.allsessions.con1(6,2) = nTrials;
lm.allsessions.con1(6,3) = nRight;
lm.allsessions.con1(6,4) = percentR;
% all other shifts
for ii = 1:5
    shift = ii*2; shiftmm = shift/10; %for naming
    name = sprintf('%d', shift);

    % creating matrix - left
    lm.allsessions.con1shift.(sprintf('left%d', shift)) = [lm.Session01.con1shift.(sprintf('left%d', shift));...
        lm.Session02.con1shift.(sprintf('left%d', shift)); lm.Session03.con1shift.(sprintf('left%d', shift));...
        lm.Session04.con1shift.(sprintf('left%d', shift))];
    % Getting key values, left shift
    nTrialsL = length(lm.allsessions.con1shift.(sprintf('left%d', shift))(:,4)); %number of trials for  the condition
    nRightL = sum(lm.allsessions.con1shift.(sprintf('left%d', shift))(:,4)==2); %number of trials they responded that the right side is longer
    percentRL = (nRightL/nTrialsL)*100;
    % creating matrix - right
    lm.allsessions.con1shift.(sprintf('right%d', shift)) = [lm.Session01.con1shift.(sprintf('right%d', shift));...
        lm.Session02.con1shift.(sprintf('right%d', shift)); lm.Session03.con1shift.(sprintf('right%d', shift));...
        lm.Session04.con1shift.(sprintf('right%d', shift))];
    % Getting key values, right shift
    nTrialsR = length(lm.allsessions.con1shift.(sprintf('right%d', shift))(:,4)); %number of trials for  the condition
    nRightR = sum(lm.allsessions.con1shift.(sprintf('right%d', shift))(:,4)==2); %number of trials they responded that the right side is longer
    percentRR = (nRightR/nTrialsR)*100;
    % Adding to matrix, left shift
    lm.allsessions.con1(k(ii),2) = nTrialsL;
    lm.allsessions.con1(k(ii),3) = nRightL;
    lm.allsessions.con1(k(ii),4) = percentRL;
    % Adding to matrix, left shift
    lm.allsessions.con1(ii+6,2) = nTrialsR;
    lm.allsessions.con1(ii+6,3) = nRightR;
    lm.allsessions.con1(ii+6,4) = percentRR;      
end

% contrast 2 - white top left
% midpoint 
lm.allsessions.con2shift.mid = [lm.Session01.con2shift.mid; lm.Session02.con2shift.mid; ...
    lm.Session03.con2shift.mid; lm.Session04.con2shift.mid];
lm.allsessions.con2(:,1) = measurements';
nTrials = length(lm.allsessions.con2shift.mid(:,4)); 
nRight = sum(lm.allsessions.con2shift.mid(:,4)==2); %number of trials they responded that the right side is longer
percentR = (nRight/nTrials)*100;
% Adding to matrix
lm.allsessions.con2(6,2) = nTrials;
lm.allsessions.con2(6,3) = nRight;
lm.allsessions.con2(6,4) = percentR;
% all other shifts
for ii = 1:5
    shift = ii*2; shiftmm = shift/10; %for naming
    name = sprintf('%d', shift);
 
    % creating matrix - left
    lm.allsessions.con2shift.(sprintf('left%d', shift)) = [lm.Session01.con2shift.(sprintf('left%d', shift));...
        lm.Session02.con2shift.(sprintf('left%d', shift)); lm.Session03.con2shift.(sprintf('left%d', shift));...
        lm.Session04.con2shift.(sprintf('left%d', shift))];
    % Getting key values, left shift
    nTrialsL = length(lm.allsessions.con2shift.(sprintf('left%d', shift))(:,4)); %number of trials for  the condition
    nRightL = sum(lm.allsessions.con2shift.(sprintf('left%d', shift))(:,4)==2); %number of trials they responded that the right side is longer
    percentRL = (nRightL/nTrialsL)*100;
     % creating matrix - right
    lm.allsessions.con2shift.(sprintf('right%d', shift)) = [lm.Session01.con2shift.(sprintf('right%d', shift));...
        lm.Session02.con2shift.(sprintf('right%d', shift)); lm.Session03.con2shift.(sprintf('right%d', shift));...
        lm.Session04.con2shift.(sprintf('right%d', shift))];
    % Getting key values, right shift
    nTrialsR = length(lm.allsessions.con2shift.(sprintf('right%d', shift))(:,4)); %number of trials for  the condition
    nRightR = sum(lm.allsessions.con2shift.(sprintf('right%d', shift))(:,4)==2); %number of trials they responded that the right side is longer
    percentRR = (nRightR/nTrialsR)*100;
    % Adding to matrix, left shift
    lm.allsessions.con2(k(ii),2) = nTrialsL;
    lm.allsessions.con2(k(ii),3) = nRightL;
    lm.allsessions.con2(k(ii),4) = percentRL;
    % Adding to matrix, left shift
    lm.allsessions.con2(ii+6,2) = nTrialsR;
    lm.allsessions.con2(ii+6,3) = nRightR;
    lm.allsessions.con2(ii+6,4) = percentRR;      
end



%% All sessions, all lines
lm.allsessions.alllines.mid = [lm.allsessions.line1.mid; lm.allsessions.line2.mid;...
     lm.allsessions.line3.mid]; %all lines
 
% All lines
lm.allsessions.res(:,1) = measurements';
nTrials = length(lm.allsessions.alllines.mid(:,4)); %number of trials for  the condition
nRight = sum(lm.allsessions.alllines.mid(:,4)==2); %number of trials they responded that the right side is longer
percentR = (nRight/nTrials)*100;
lm.allsessions.res(6,2) = nTrials;
lm.allsessions.res(6,3) = nRight;
lm.allsessions.res(6,4) = percentR;

% All other offsets
for ii = 1:5 %number of offsets for left/right
    shift = ii*2; shiftmm = shift/10; %for naming
    name = sprintf('%d', shift);
    % All lines
    % Matrices
    lm.allsessions.alllines.(sprintf('left%d', shift)) = [lm.allsessions.line1.(sprintf('left%d', shift));...
        lm.allsessions.line2.(sprintf('left%d', shift)); lm.allsessions.line3.(sprintf('left%d', shift))]; %left
    lm.allsessions.alllines.(sprintf('right%d', shift)) = [lm.allsessions.line1.(sprintf('right%d', shift));...
        lm.allsessions.line2.(sprintf('right%d', shift)); lm.allsessions.line3.(sprintf('right%d', shift))]; %right
    % Key values, left
    nTrialsL = length(lm.allsessions.alllines.(sprintf('left%d', shift))(:,4)); %number of trials for  the condition
    nRightL = sum(lm.allsessions.alllines.(sprintf('left%d', shift))(:,4)==2); %number of trials they responded that the right side is longer
    percentRL = (nRightL/nTrialsL)*100;
    % Key values, right
    nTrialsR = length(lm.allsessions.alllines.(sprintf('right%d', shift))(:,4)); %number of trials for  the condition
    nRightR = sum(lm.allsessions.alllines.(sprintf('right%d', shift))(:,4)==2); %number of trials they responded that the right side is longer
    percentRR = (nRightR/nTrialsR)*100;
    % Matrix left
    lm.allsessions.res(k(ii),2) = nTrialsL;
    lm.allsessions.res(k(ii),3) = nRightL;
    lm.allsessions.res(k(ii),4) = percentRL;
    % Matrix right
    lm.allsessions.res(ii+6,2) = nTrialsR;
    lm.allsessions.res(ii+6,3) = nRightR;
    lm.allsessions.res(ii+6,4) = percentRR;
end

    %% Lapse rate for LM task 
    % calculating response error for the lines that are clearly longer on l/r
    % (10mm shift)
    errorL = (sum(lm.allsessions.alllines.left10(:,4)==2)/...
        length(lm.allsessions.alllines.left10(:,4)))*100; %left 10mm shift, finding incorrect responses
    errorR = (sum(lm.allsessions.alllines.right10(:,4)==1)/...
        length(lm.allsessions.alllines.right10(:,4)))*100; %right 10mm shift, finding incorrect responses
    errors = [errorL, errorR];
    lm.lapse = mean(errors); %lapse rate calculated in percentage

    %% Plotting LM task
    cd(dirVis); %analysis directory for saving figures
    % Plot percent left side longer for each session
    % For each line 
    asym = lm.Session01.line1.res(:,1); %for labelling
    for i = 1:3 %number of lines
        line = sprintf('line%d', i);
        figure(i+5) 
        plot(asym, lm.Session01.(sprintf('%s', line)).res(:,4));
        hold on
        plot(asym, lm.Session02.(sprintf('%s', line)).res(:,4));
        hold on
        plot(asym, lm.Session03.(sprintf('%s', line)).res(:,4));
        hold on
        plot(asym, lm.Session04.(sprintf('%s', line)).res(:,4));
        % Adding the shifts as x-axis tick labels
        ax = gca;
        set(ax, 'Xtick', asym);
        xlim([min(asym) max(asym)]); ylim([0 100]);
        % Adding horizontal line at y = 50
        ymid = max(ylim)/2;
        hold on
        plot(xlim, [1,1]*ymid, '--k')
        % Naming
        legend('Sess1', 'Sess2', 'Sess3', 'Sess4');
        ylabel('Left-side perceived as longer (%)');
        xlabel('Stimulus asymmetry (mm)');
        saveas(figure(i+4), sprintf('%s_LM%s_raw.jpg', ppID, line));
    end

    % For all lines
    figure(9)
    plot(asym, lm.Session01.res(:,4));
    hold on
    plot(asym, lm.Session02.res(:,4));
    hold on
    plot(asym, lm.Session03.res(:,4));
    hold on
    plot(asym, lm.Session04.res(:,4));
    % Adding the shifts as x-axis tick labels
    ax = gca;
    set(ax, 'Xtick', asym);
    xlim([min(asym) max(asym)]); ylim([0 100]);
    % Adding horizontal line at y = 50
    ymid = max(ylim)/2;
    hold on
    plot(xlim, [1,1]*ymid, '--k')
    % Naming
    legend('Sess1', 'Sess2', 'Sess3', 'Sess4');
    ylabel('Left-side perceived as longer (%)');
    xlabel('Stimulus asymmetry (mm)');
    saveas(figure(9), sprintf('%s_LMsess_raw.jpg', ppID));

    % For all sessions
    figure(10)
    plot(asym, lm.allsessions.line1.res(:,4));
    hold on
    plot(asym, lm.allsessions.line2.res(:,4));
    hold on
    plot(asym, lm.allsessions.line3.res(:,4));
    % Adding the shifts as x-axis tick labels
    ax = gca;
    set(ax, 'Xtick', asym);
    xlim([min(asym) max(asym)]); ylim([0 100]);
    % Adding horizontal line at y = 50
    ymid = max(ylim)/2;
    hold on
    plot(xlim, [1,1]*ymid, '--k')
    % Naming
    legend('10cm line', '20cm line', '30cm line');
    ylabel('Left-side perceived as longer (%)');
    xlabel('Stimulus asymmetry (mm)');
    saveas(figure(10), sprintf('%s_LMlines_raw.jpg', ppID));

    %% Adding psychometric functions to here
    % Use the demo and have a look at Abi's code

    %% Analyse LM2 data
    % Take percentage top-line longer responses for each shift in mm (of the
    % top line)
    % Group ing into line length
    for i = 1:length(nSessions)
        session = sprintf('Session%0*d',2,nSessions(i));
        cd(dirSess); %directing to current session folder
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

       %% Find percentage likelihood left shifted line is perceived as longer
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
        lm2.(sprintf('%s', session)).line1.per(6, 2) = perMid;
        perMid = sum((lm2.(sprintf('%s', session)).line2.mid(:,9))/...
            length(lm2.(sprintf('%s', session)).line2.mid(:,9)))*100; %line 2
        lm2.(sprintf('%s', session)).line2.per(6, 2) = perMid;
        perMid = sum((lm2.(sprintf('%s', session)).line3.mid(:,9))/...
            length(lm2.(sprintf('%s', session)).line3.mid(:,9)))*100; %line 3
        lm2.(sprintf('%s', session)).line3.per(6, 2) = perMid;

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
                lm2.(sprintf('%s', session)).line1.(sprintf('left%d', shift))(j,8)); %left, is equal
                lm2.(sprintf('%s', session)).line1.(sprintf('right%d', shift))(j,9) = ...
                    ne(lm2.(sprintf('%s', session)).line1.(sprintf('right%d', shift))(j,4),...
                lm2.(sprintf('%s', session)).line1.(sprintf('right%d', shift))(j,8)); %right, is NOT equal (don't expect response to match RHS shifted line) 
                % Line 2
                lm2.(sprintf('%s', session)).line2.(sprintf('left%d', shift))(j,9) = ...
                    isequal(lm2.(sprintf('%s', session)).line2.(sprintf('left%d', shift))(j,4),...
                lm2.(sprintf('%s', session)).line2.(sprintf('left%d', shift))(j,8)); %left, is equal
                lm2.(sprintf('%s', session)).line2.(sprintf('right%d', shift))(j,9) = ...
                    isequal(lm2.(sprintf('%s', session)).line2.(sprintf('right%d', shift))(j,4),...
                lm2.(sprintf('%s', session)).line2.(sprintf('right%d', shift))(j,8)); %right, is NOT equal 
                % Line 3
                lm2.(sprintf('%s', session)).line3.(sprintf('left%d', shift))(j,9) = ...
                    isequal(lm2.(sprintf('%s', session)).line3.(sprintf('left%d', shift))(j,4),...
                lm2.(sprintf('%s', session)).line3.(sprintf('left%d', shift))(j,8)); %left, is equal
                lm2.(sprintf('%s', session)).line3.(sprintf('right%d', shift))(j,9) = ...
                    ne(lm2.(sprintf('%s', session)).line3.(sprintf('right%d', shift))(j,4),...
                lm2.(sprintf('%s', session)).line3.(sprintf('right%d', shift))(j,8)); %right, is not equal
            end
        end

    end
    %% Plotting LM2 task
    %% Lapse rate for LM2 task
    %% Close and save
    close all
    cd(dirVis)
    save(matfilename, 'mlb', 'lm', 'lm2');
end