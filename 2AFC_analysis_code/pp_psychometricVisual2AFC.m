%% AG. Mitchell - 23.01.19
%% Analysis for landmarks 2AFC task
% Takes the data from the landmarks 2AFC task and calculates percentage of
% 'top line longer' responses for each left/right shift.
% Calculates psychometric curve (cumulative normal) based on this data

%% Load data
%% Edited demo to fit bias study data (LM task) AG.Mitchell 22.11.18
clear all;      %Clear all existing variables from memory

tic
% Load in participant data
nParticipants = [1:19,21:24,26:30];
%nParticipants = 3; %for testing
for p = 1:length(nParticipants)  
    ppID = sprintf('P%0*d',2,nParticipants(p));
    %ppID = input('Participant ID? ', 's'); %for use when navigating files
    visfilename = sprintf('%s_visualanalysisStart.mat', ppID);
    matfilename = sprintf('%s_visualanalysis2AFC.mat', ppID);
    nSessions = 1:4; %vector number of sessions each participant does
    % Directory
    dirBias = ('M:\Alex_Files\Experiments\Bias\'); %subject to change depending on where you analyse
    dirPP = [dirBias 'Data' filesep ppID]; %participant directory
    % Navigating to analysis folder
    cd(dirPP)
    dirAna = [dirPP filesep 'Analysis' filesep];
    dirVis = [dirAna 'Visual' filesep];
    cd(dirVis)
    load(visfilename)

%% Calculate proportion response
% Calculating the proportion of responses that said the line shifted to the
% RHS was longer
% Use the 'all shift' for this
    for i = 1:length(nSessions)
        session = sprintf('Session%0*d',2,nSessions(i));
        dirSess = [dirPP filesep session];
        cd(dirSess); %directing to current session folder
        % Getting measurements of shift
        measurements = -10:2:10;
        lm2.(sprintf('%s', session)).allshift.per(:,1) = measurements';
        
        % Calculating vectors of equal values between 'line to shift left' and
        % 'response' - used to be able to calculate percentage of
        % 'left-shitfted line is longer' responses
        lm2mat = lm2.(sprintf('%s', session)).matrix;
        lm2.(sprintf('%s', session)).allshift.midL = lm2mat(find(lm2mat(:,2)== 0),:); %midpoint
        lm2.(sprintf('%s', session)).allshift.midR = lm2mat(find(lm2mat(:,3 )== 0),:);
        % Midpoint lines - for each line
        % Left mid
        for j = 1:length(lm2.(sprintf('%s', session)).allshift.midL(:,4))
            lm2.(sprintf('%s', session)).allshift.midL(j,9) = ...
                ~isequal(lm2.(sprintf('%s', session)).allshift.midL(j,4), lm2.(sprintf('%s', session)).allshift.midL(j,8));
        end
        % Right mid
        for j = 1:length(lm2.(sprintf('%s', session)).allshift.midR(:,4))
            lm2.(sprintf('%s', session)).allshift.midR(j,9) = ...
                ~isequal(lm2.(sprintf('%s', session)).allshift.midR(j,4), lm2.(sprintf('%s', session)).allshift.midR(j,8));
        end
        % Calculating the percentage for the midpoint
        perMidL = sum((lm2.(sprintf('%s', session)).allshift.midL(:,9))/...
            length(lm2.(sprintf('%s', session)).allshift.midL(:,9)))*100;
        perMidR = sum((lm2.(sprintf('%s', session)).allshift.midR(:,9))/...
            length(lm2.(sprintf('%s', session)).allshift.midR(:,9)))*100;
        perMid = (perMidL+perMidR)/2;
        
        lm2.(sprintf('%s', session)).allshift.per(6,2) = perMid;
        
        % Do the same for all the shifts now
        k = 1:5; k = flipud(k');
        for ii = 1:5 %number of shifts
            shift = ii*2; shiftmm = shift/10; %for naming
            name = sprintf('%d', shift);
            for j = 1:length(lm2.(sprintf('%s', session)).allshift.(sprintf('left%d', shift))(:,4))
                lm2.(sprintf('%s', session)).allshift.(sprintf('left%d', shift))(j,9) = ...
                        ~isequal(lm2.(sprintf('%s', session)).allshift.(sprintf('left%d', shift))(j,4),...
                lm2.(sprintf('%s', session)).allshift.(sprintf('left%d', shift))(j,8)); 
                lm2.(sprintf('%s', session)).allshift.(sprintf('right%d', shift))(j,9) = ...
                    ~isequal(lm2.(sprintf('%s', session)).allshift.(sprintf('right%d', shift))(j,4),...
                lm2.(sprintf('%s', session)).allshift.(sprintf('right%d', shift))(j,8)); 
            end
            perLeft = sum((lm2.(sprintf('%s', session)).allshift.(sprintf('left%d', shift))(:,9))/...
            length(lm2.(sprintf('%s', session)).allshift.(sprintf('left%d', shift))(:,9)))*100;
            perRight = sum((lm2.(sprintf('%s', session)).allshift.(sprintf('right%d', shift))(:,9))/...
            length(lm2.(sprintf('%s', session)).allshift.(sprintf('right%d', shift))(:,9)))*100;
        
            lm2.(sprintf('%s', session)).allshift.per(k(ii),2) = perLeft;
            lm2.(sprintf('%s', session)).allshift.per((ii+6),2) = perRight;
        end
        
        lm2psych.(sprintf('%s', session)) = lm2.(sprintf('%s', session)).allshift.per;
    end
    
    %average across all sessions
    allSessions = [lm2psych.Session01(:,2), lm2psych.Session02(:,2), lm2psych.Session03(:,2),...
        lm2psych.Session04(:,2)];
    asym = lm2psych.Session01(:,1); %asymmtry for the x-axis
    lm2psych.allSessions = [asym, mean(allSessions,2)]; % mean all sessions
    
    %Saving
    cd(dirVis); save(matfilename, 'lm2psych');
    
    %% Creating a plot
    % All sessions
    asym = [-10:2:0]; 
    % All sessions
    fig1name = sprintf('%slm2raw_Sessions', ppID);
    pdfFileName = strcat(fig1name, '.pdf');
    for i = 1:length(nSessions)
        session = sprintf('Session%0*d',2,nSessions(i));
        sess = lm2psych.(sprintf('%s', session))(1:6,2);
        figure(1)
        plot(asym, sess, 'LineWidth', 1);
        hold on
    end 
    % Adding the shifts as x-axis tick labels
    figure(1)
    ax = gca;
    set(ax, 'Xtick', asym);
    xlim([min(asym) max(asym)]); ylim([0 100]);
    % Adding horizontal line at y = 50
    ymid = max(ylim)/2;
    hold on
    plot(xlim, [1,1]*ymid, '--k')
    % Naming
    legend('Sess1', 'Sess2', 'Sess3', 'Sess4');
    ylabel('Perceneage right shifted line perceived as longer');
    xlabel('Stimulus asymmetry (mm)'); title('Visual');
    saveas(gcf, pdfFileName);
    
    % Average session
    fig2name = sprintf('%slm2raw_meanSessions', ppID);
    pdfFileName = strcat(fig2name, '.pdf');
    
    figure(2)
    plot(asym, lm2psych.allSessions(1:6,2), 'LineWidth', 1.5);
    % Adding the shifts as x-axis tick labels
    ax = gca;
    set(ax, 'Xtick', asym);
    xlim([min(asym) max(asym)]); ylim([0 100]);
    % Adding horizontal line at y = 50
    ymid = max(ylim)/2;
    hold on
    plot(xlim, [1,1]*ymid, '--k')
    ylabel('Perceneage right shifted line perceived as longer');
    xlabel('Stimulus asymmetry (mm)'); title('Visual');
    saveas(gcf, pdfFileName);
   
%% Psychometric curve fitting
close all
end