%
%PAL_PFML_Demo  Demonstrates use of Palamedes functions to (1) fit a
%Psychometric Function to some data using a Maximum Likelihood criterion, 
%(2) determine standard errors of free parameters using a bootstrap
%procedure and (3) determine the goodness-of-fit of the fit.
%
%Demonstrates basic usage of Palamedes functions:
%-PAL_PFML_Fit
%-PAL_PFML_BootstrapParametric
%-PAL_PFML_BootstrapNonParametric
%-PAL_PFML_GoodnessOfFit
%secondary:
%-PAL_Logistic
%
%More information on any of these functions may be found by typing
%help followed by the name of the function. e.g., help PAL_PFML_Fit
%
%NP (08/23/2009)

%% Edited demo to fit bias study data (LM task) AG.Mitchell 22.11.18
clear all;      %Clear all existing variables from memory

tic
% Load in participant data
nParticipants = [1:17,21,22];
for p = 1:length(nParticipants)  
    ppID = sprintf('P%0*d',2,nParticipants(p));
    %ppID = input('Participant ID? ', 's'); %for use when navigating files
    visfilename = sprintf('%s_visualanalysisStart.mat', ppID);
    matfilename = sprintf('%s_visualanalysis.mat', ppID);
    nSessions = 1:4; %vector number of sessions each participant does
    % Directory
    dirBias = ('C:\Users\Experimenter\Documents\Experiments2018\Bias'); %subject to change depending on where you analyse
    dirPP = [dirBias filesep ppID]; %participant directory
    % Making new anaysis folder for saving
    cd(dirPP)
    mkdir Analysis;
    dirAna = [dirPP filesep 'Analysis' filesep];
    cd(dirAna)
    mkdir Visual
    dirVis = [dirAna 'Visual' filesep];
    cd(dirVis)
    load(visfilename)

    %% Variables for curve fitting
    ParOrNonPar = 2; %non-parametric design
    nSims = 1000; %number of sims for all bootstraps and goodness of fit
    lm.lapseP = lm.lapse/100; %proportion correct, not percentage
    paramsFree = [1 1 0 0];  %1: free parameter, 0: fixed parameter

    %Use the Cumulative normal function
    PF = @PAL_CumulativeNormal;  %Alternatives: PAL_Gumbel, PAL_Weibull,
                         %PAL_Quick, PAL_logQuick,
                         %PAL_CumulativeNormal, PAL_HyperbolicSecant


    fitLapseRate = false;

    %For this set of stim values the guess rate is equal to the lapse rate.  
    %This can be kept true for this data. 
    gammaEqLambda = true;

    %This is just to show how it would be done.
    %Really for lapse estimation to work well should do it with a multiple fit procedure.
    if fitLapseRate %Try and fit the lapse rate. 
        %Fit lapse rate:
        searchGrid.gamma = [0:.005:.1];     %type help PAL_PFML_Fit for more information
        searchGrid.lambda = [0:.005:.1];
        paramsFree = [1 1 1 1];  %1: free parameter, 0: fixed parameter
    else %Don't fit lapse rate and fix it to a value
        searchGrid.gamma = lm.lapseP;  %This is the lapse rate for negative stim %values
        searchGrid.lambda = lm.lapseP;  %this is lapse rate for the participant, calculated here
        paramsFree = [1 1 0 0];
    end

    %% Analysis per session
    for i = 1:length(nSessions)
        session = sprintf('Session%0*d',2,nSessions(i));
        StimLevels = lm.(sprintf('%s', session)).res(:,1)'; %stimulus asymmetry
        % Number of 'right-side longer' responses for each stim level
        NumPos = lm.(sprintf('%s', session)).res(:,3)';             
        % Number of trials at each entry of 'StimLevels'
        % Percentage so out of 100
        OutOfNum = lm.(sprintf('%s', session)).res(:,2)';

        %Parameter grid defining parameter space through which to perform a
        %brute-force search for values to be used as initial guesses in iterative
        %parameter search.
        searchGrid.alpha = linspace(min(StimLevels), max(StimLevels), 101);
        searchGrid.beta = linspace(0,30/max(StimLevels),101); %slope

        % Perform fit - Cumulative Normal
        disp('Fitting function.....');
        [paramsValues LL exitflag] = PAL_PFML_Fit(StimLevels,NumPos, ...
        OutOfNum,searchGrid,paramsFree,PF,...
        'lapseLimits',[0 1],'gammaEQlambda',gammaEqLambda)

        % Getting standard error
        if ParOrNonPar == 1
            [SD paramsSim LLSim converged] = PAL_PFML_BootstrapParametric(...
                StimLevels, OutOfNum, paramsValues, paramsFree, nSims, PF, ...
                'searchGrid', searchGrid);
        else
            [SD paramsSim LLSim converged] = PAL_PFML_BootstrapNonParametric(...
                StimLevels, NumPos, OutOfNum, [], paramsFree, nSims, PF,...
                'searchGrid',searchGrid);
        end

        [Dev pDev] = PAL_PFML_GoodnessOfFit(StimLevels, NumPos, OutOfNum, ...
        paramsValues, paramsFree, nSims, PF, 'searchGrid', searchGrid);

        %Create simple plot
        ProportionCorrectObserved=NumPos./OutOfNum; 
        StimLevelsFineGrain=[min(StimLevels):max(StimLevels)./1000:max(StimLevels)];
        ProportionCorrectModel = PF(paramsValues,StimLevelsFineGrain);

        cd(dirVis)
        figure();
        axes
        hold on
        plot(StimLevelsFineGrain,ProportionCorrectModel,'-','color',[.6 0 .2],'linewidth',2);
        plot(StimLevels,ProportionCorrectObserved,'k.','markersize',20);
        %observedHandle = plot2afc(StimLevels,NumPos,OutOfNum);
        %set(observedHandle,'color','k','linewidth',1.5)
        set(gca, 'fontsize',14);
        set(gca, 'Xtick',StimLevels);
        axis([min(StimLevels) max(StimLevels) 0 1]);
        xlabel('Stimulus Intensity');
        ylabel('Proportion right-side longer');
        figFileName = strcat(ppID, '_', 'pFit', session, '.pdf');
        saveas(gcf, figFileName);

        % Getting the information that we need at 50% 
        stim50right = PAL_CumulativeNormal(paramsValues, 0.5, 'Inverse');
        slope50thresh = PAL_CumulativeNormal(paramsValues, stim50right, 'Derivative');

        for iBoot = 1:nSims
            boot50thresh(iBoot) = PAL_CumulativeNormal(paramsSim(iBoot,:), 0.5, 'Inverse');
            bootSlope50thresh(iBoot)= PAL_CumulativeNormal(paramsSim(iBoot,:), boot50thresh(iBoot), 'Derivative');
        end
        % Getting the standard error for the psychometric data
        thresholdSE = std(boot50thresh);
        slopeSE = std(bootSlope50thresh);

        % Getting confidence intervals for the function
        sortedThresholdSim = sort(boot50thresh);
        sortedSlopeSim = sort(bootSlope50thresh);
        thresholdCI = [sortedThresholdSim(25) sortedThresholdSim(nSims-25)];
        slopeCI = [sortedSlopeSim(25) sortedSlopeSim(nSims-25)];

        lm.pfits.(sprintf('%s', session)).stim50right = stim50right;
        lm.pfits.(sprintf('%s', session)).slope50thresh = slope50thresh;
        lm.pfits.(sprintf('%s', session)).threshCI = thresholdCI;
        lm.pfits.(sprintf('%s', session)).slopeCI = slopeCI;
        lm.pfits.(sprintf('%s', session)).params.values = paramsValues;
        lm.pfits.(sprintf('%s', session)).params.SE = SD;
        lm.pfits.(sprintf('%s', session)).params.Dev = Dev;
        lm.pfits.(sprintf('%s', session)).params.pDev = pDev;
        lm.pfits.(sprintf('%s', session)).model = ProportionCorrectModel; %last one for plotting :)
    end

    % Plot of all session fits
    figure();
    axes
    hold on
    plot(StimLevelsFineGrain,lm.pfits.Session01.model,'-','color',[.6 0 .1],'linewidth',2);
    hold on
    plot(StimLevelsFineGrain,lm.pfits.Session02.model,'-','color',[0 .6 .1],'linewidth',2);
    hold on
    plot(StimLevelsFineGrain,lm.pfits.Session03.model,'-','color',[0 .1 .6],'linewidth',2);
    hold on
    plot(StimLevelsFineGrain,lm.pfits.Session04.model,'-','color',[.2 .2 .2],'linewidth',2);
    set(gca, 'fontsize',14);
    set(gca, 'Xtick',StimLevels);
    axis([min(StimLevels) max(StimLevels) 0 1]);
    xlabel('Stimulus Intensity');
    ylabel('Proportion right-side longer');
    legend(' 1', ' 2', ' 3', ' 4', 'Position', [350 90 0.2 0.1]);
    figFileName = strcat(ppID, '_', 'pFit_allsess', '.pdf');
    saveas(gcf, figFileName);

    %% All sessions and all lines
    % Getting correct data
    StimLevels = lm.allsessions.res(:,1)'; %stimulus asymmetry
    NumPos = lm.allsessions.res(:,3)';             
    OutOfNum = lm.allsessions.res(:,2)';

    % Search grid
    searchGrid.alpha = linspace(min(StimLevels), max(StimLevels), 101);
    searchGrid.beta = linspace(0,30/max(StimLevels),101); %slope

    % Perform fit - Cumulative Normal
    disp('Fitting function.....');
    [paramsValues LL exitflag] = PAL_PFML_Fit(StimLevels,NumPos, ...
        OutOfNum,searchGrid,paramsFree,PF,...
        'lapseLimits',[0 1],'gammaEQlambda',gammaEqLambda)

    % Getting standard error
    if ParOrNonPar == 1
        [SD paramsSim LLSim converged] = PAL_PFML_BootstrapParametric(...
            StimLevels, OutOfNum, paramsValues, paramsFree, nSims, PF, ...
            'searchGrid', searchGrid);
    else
        [SD paramsSim LLSim converged] = PAL_PFML_BootstrapNonParametric(...
            StimLevels, NumPos, OutOfNum, [], paramsFree, nSims, PF,...
            'searchGrid',searchGrid);
    end


    disp('Determining Goodness-of-fit.....');

    [Dev pDev] = PAL_PFML_GoodnessOfFit(StimLevels, NumPos, OutOfNum, ...
        paramsValues, paramsFree, nSims, PF, 'searchGrid', searchGrid);

    %Create simple plot
    ProportionCorrectObserved=NumPos./OutOfNum; 
    StimLevelsFineGrain=[min(StimLevels):max(StimLevels)./1000:max(StimLevels)];
    ProportionCorrectModel = PF(paramsValues,StimLevelsFineGrain);

    cd(dirVis)
    figure();
    axes
    hold on
    plot(StimLevelsFineGrain,ProportionCorrectModel,'-','color',[0 .6 .2],'linewidth',2);
    plot(StimLevels,ProportionCorrectObserved,'k.','markersize',20);
    %observedHandle = plot2afc(StimLevels,NumPos,OutOfNum);
    %set(observedHandle,'color','k','linewidth',1.5)
    set(gca, 'fontsize',14);
    set(gca, 'Xtick',StimLevels);
    axis([min(StimLevels) max(StimLevels) 0 1]);
    xlabel('Stimulus Intensity');
    ylabel('Proportion right-side longer');
    figFileName = strcat(ppID, '_', 'pFitAll',  '.pdf');
    saveas(gcf, figFileName);

    % Getting the information that we need at 50% 
    stim50right = PAL_CumulativeNormal(paramsValues, 0.5, 'Inverse');
    slope50thresh = PAL_CumulativeNormal(paramsValues, stim50right, 'Derivative');

    for iBoot = 1:nSims
        boot50thresh(iBoot) = PAL_CumulativeNormal(paramsSim(iBoot,:), 0.5, 'Inverse');
        bootSlope50thresh(iBoot)= PAL_CumulativeNormal(paramsSim(iBoot,:), boot50thresh(iBoot), 'Derivative');
    end
    % Getting the standard error for the psychometric data
    thresholdSE = std(boot50thresh);
    slopeSE = std(bootSlope50thresh);

    % Getting confidence intervals for the function
    sortedThresholdSim = sort(boot50thresh);
    sortedSlopeSim = sort(bootSlope50thresh);
    thresholdCI = [sortedThresholdSim(25) sortedThresholdSim(nSims-25)];
    slopeCI = [sortedSlopeSim(25) sortedSlopeSim(nSims-25)];

    % Adding these functions to the data matrix
        % Adding these functions to the data matrix
    lm.pfits.all.stim50right = stim50right;
    lm.pfits.all.slope50thresh = slope50thresh;
    lm.pfits.all.threshCI = thresholdCI;
    lm.pfits.all.slopeCI = slopeCI;
    lm.pfits.all.params.values = paramsValues;
    lm.pfits.all.params.SE = SD;
    lm.pfits.all.params.Dev = Dev;
    lm.pfits.all.params.pDev = pDev;
    lm.pfits.all.model = ProportionCorrectModel;

    close all
    cd(dirVis)
    save(matfilename, 'lm', 'lm2', 'mlb');

end
toc