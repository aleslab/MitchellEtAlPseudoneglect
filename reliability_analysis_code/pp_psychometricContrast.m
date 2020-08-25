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
filePath = cd;
[dirBias, name, ext] = fileparts(filePath); %subject to change depending on where you analyse
dirData = [dirBias filesep 'Data'];
    
nParticipants = [1:19,21:24,26:30];
%nParticipants = 11; %for testing
for p = 1:length(nParticipants)  
    ppID = sprintf('P%0*d',2,nParticipants(p));
    %ppID = input('Participant ID? ', 's'); %for use when navigating files
    visfilename = sprintf('%s_visualanalysisStart.mat', ppID);
    matfilename = sprintf('%s_visualanalysisContrast.mat', ppID);
    nSessions = 1:2; %vector number of sessions each participant does
    % Directory
    dirPP = [dirData filesep ppID]; %participant directory
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
    paramsFree = [1 1 1 1];  %1: free parameter, 0: fixed parameter

    %Use the Cumulative normal function
    PF = @PAL_CumulativeNormal;  %Alternatives: PAL_Gumbel, PAL_Weibull,
                         %PAL_Quick, PAL_logQuick,
                         %PAL_CumulativeNormal, PAL_HyperbolicSecant


    fitLapseRate = true;

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
        con = sprintf('con%d',nSessions(i));
        StimLevels = lm.allsessions.(sprintf('%s', con))(:,1)'; %stimulus asymmetry
        % Number of 'right-side longer' responses for each stim level
        NumPos = lm.allsessions.(sprintf('%s', con))(:,3)';             
        % Number of trials at each entry of 'StimLevels'
        % Percentage so out of 100
        OutOfNum = lm.allsessions.(sprintf('%s', con))(:,2)';
        options = PAL_minimize('options'); %for lapse rate fitting

        %Parameter grid defining parameter space through which to perform a
        %brute-force search for values to be used as initial guesses in iterative
        %parameter search.
        searchGrid.alpha = linspace(min(StimLevels), max(StimLevels), 101);
        searchGrid.beta = linspace(0,30/max(StimLevels),101); %slope

        % Perform fit - Cumulative Normal
        disp('Fitting function.....');
        [paramsValues LL exitflag] = PAL_PFML_Fit(StimLevels,NumPos, ...
        OutOfNum,searchGrid,paramsFree,PF,...
        'lapseLimits',[0 1],'gammaEQlambda', gammaEqLambda)

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
        paramsValues, paramsFree, nSims, PF, 'lapseLimits', [0 1], 'gammaEQlambda', gammaEqLambda,...
        'searchGrid', searchGrid);

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
        set(gca, 'fontsize',16);
        set(gca, 'Xtick',StimLevels);
        axis([min(StimLevels) max(StimLevels) 0 1]);
        xlabel('Stimulus Asymmetry (mm)');
        ylabel('Proportion right-side longer');
        figFileName = strcat(ppID, '_', 'pFit', con, '.pdf');
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

        lm.pfits.(sprintf('%s', con)).stim50right = stim50right;
        lm.pfits.(sprintf('%s', con)).slope50thresh = slope50thresh;
        lm.pfits.(sprintf('%s', con)).threshCI = thresholdCI;
        lm.pfits.(sprintf('%s', con)).slopeCI = slopeCI;
        lm.pfits.(sprintf('%s', con)).params.values = paramsValues;
        lm.pfits.(sprintf('%s', con)).params.SE = SD;
        lm.pfits.(sprintf('%s', con)).params.Dev = Dev;
        lm.pfits.(sprintf('%s', con)).params.pDev = pDev;
        lm.pfits.(sprintf('%s', con)).model = ProportionCorrectModel; %last one for plotting :)
    end

    % Plot of all session fits
    figure();
    axes
    hold on
    plot(StimLevelsFineGrain,lm.pfits.con1.model,'-','color',[.6 0 .1],'linewidth',2);
    hold on
    plot(StimLevelsFineGrain,lm.pfits.con2.model,'-','color',[0 .6 .1],'linewidth',2);
    set(gca, 'fontsize',15);
    set(gca, 'Xtick',StimLevels);
    axis([min(StimLevels) max(StimLevels) 0 1]);
    xlabel('Stimulus Asymmetry (mm)');
    ylabel('Proportion right-side longer');
    legend('con1', 'con2', 'Position', [350 90 0.2 0.1]);
    figFileName = strcat(ppID, '_', 'pFit_allsess', '.pdf');
    saveas(gcf, figFileName);

   
    close all
    cd(dirVis)
    save(matfilename, 'lm', 'mlb');

end
toc