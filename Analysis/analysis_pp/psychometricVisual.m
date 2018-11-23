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
ppID = input('Participant ID? ', 's'); %for use when navigating files
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
load(matfilename)

%% Variables for curve fitting
ParOrNonPar = 2; %non-parametric design
nSims = 1000; %number of sims for all bootstraps and goodness of fit
lm.lapse = lm.lapse/100; %proportion correct, not percentage
paramsFree = [1 1 0 0];  %1: free parameter, 0: fixed parameter

%Use the Cumulative normal function
PF = @PAL_CumulativeNormal;  %Alternatives: PAL_Gumbel, PAL_Weibull,
                     %PAL_Quick, PAL_logQuick,
                     %PAL_CumulativeNormal, PAL_HyperbolicSecant

%% All sessions and all lines
% Stimulus intensities, this is the 'x-axis' of the plot (aka. IV)
StimLevels = lm.allsessions.res(:,1)'; %stimulus asymmetry
% Number of 'right-side longer' responses for each stim level
NumPos = lm.allsessions.res(:,3)';             
% Number of trials at each entry of 'StimLevels'
% Percentage so out of 100
OutOfNum = lm.allsessions.res(:,2)';

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
    searchGrid.gamma = lm.lapse;  %This is the lapse rate for negative stim %values
    searchGrid.lambda = lm.lapse;  %this is lapse rate for the participant, calculated here
    paramsFree = [1 1 0 0];
end

%Parameter grid defining parameter space through which to perform a
%brute-force search for values to be used as initial guesses in iterative
%parameter search.
searchGrid.alpha = linspace(min(StimLevels), max(StimLevels), 101);
searchGrid.beta = linspace(0,30/max(StimLevels),101); %slope

%% Perform fit - Cumulative Normal
disp('Fitting function.....');
[paramsValues LL exitflag] = PAL_PFML_Fit(StimLevels,NumPos, ...
    OutOfNum,searchGrid,paramsFree,PF,...
    'lapseLimits',[0 1],'gammaEQlambda',gammaEqLambda)

disp('done:')
message = sprintf('Threshold estimate: %6.4f',paramsValues(1));
disp(message);
message = sprintf('Slope estimate: %6.4f\r',paramsValues(2));
disp(message);

disp('Determining standard errors.....');        

if ParOrNonPar == 1
    [SD paramsSim LLSim converged] = PAL_PFML_BootstrapParametric(...
        StimLevels, OutOfNum, paramsValues, paramsFree, nSims, PF, ...
        'searchGrid', searchGrid);
else
    [SD paramsSim LLSim converged] = PAL_PFML_BootstrapNonParametric(...
        StimLevels, NumPos, OutOfNum, [], paramsFree, nSims, PF,...
        'searchGrid',searchGrid);
end

disp('done:');
message = sprintf('Standard error of Threshold: %6.4f',SD(1));
disp(message);
message = sprintf('Standard error of Slope: %6.4f\r',SD(2));
disp(message);


disp('done:');
message = sprintf('Standard error of Threshold: %6.4f',SD(1));
disp(message);
message = sprintf('Standard error of Slope: %6.4f\r',SD(2));
disp(message);

%Distribution of estimated slope parameters for simulations will be skewed
%(type: hist(paramsSim(:,2),40) to see this). However, distribution of
%log-transformed slope estimates will be approximately symmetric
%[type: hist(log10(paramsSim(:,2),40)]. This might motivate using 
%log-scale for slope values (uncomment next three lines to put on screen):
% SElog10slope = std(log10(paramsSim(:,2)));
% message = ['Estimate for log10(slope): ' num2str(log10(paramsValues(2))) ' +/- ' num2str(SElog10slope)];
% disp(message);

%Number of simulations to perform to determine Goodness-of-Fit
nSims=1000;

disp('Determining Goodness-of-fit.....');

[Dev pDev] = PAL_PFML_GoodnessOfFit(StimLevels, NumPos, OutOfNum, ...
    paramsValues, paramsFree, nSims, PF, 'searchGrid', searchGrid);

disp('done:');

%Put summary of results on screen
message = sprintf('Deviance: %6.4f',Dev);
disp(message);
message = sprintf('p-value: %6.4f',pDev);
disp(message);
 
%Create simple plot
ProportionCorrectObserved=NumPos./OutOfNum; 
StimLevelsFineGrain=[min(StimLevels):max(StimLevels)./1000:max(StimLevels)];
ProportionCorrectModel = PF(paramsValues,StimLevelsFineGrain);
 
figure('name','Maximum Likelihood Psychometric Function Fitting');
axes
hold on
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','color',[0 .7 .2],'linewidth',2);
 
%plot(StimLevels,ProportionCorrectObserved,'k.','markersize',40);
observedHandle = plot2afc(StimLevels,NumPos,OutOfNum);
set(observedHandle,'color','k','linewidth',1.5)

set(gca, 'fontsize',14);
set(gca, 'Xtick',StimLevels);
axis([min(StimLevels) max(StimLevels) 0 1]);
xlabel('Stimulus Intensity');
ylabel('Proportion right-side longer');

% Getting the information that we need at 50% 
stim50Right = PAL_CumulativeNormal(paramsValues, 0.5, 'Inverse');
slope50thresh = PAL_CumulativeNormal(paramsValues, stim50Right, 'Derivative');

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

% Adding these functions to the data

toc