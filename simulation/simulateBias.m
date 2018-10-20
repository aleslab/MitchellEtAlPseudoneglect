function results=simulateBias()
%To simulate the data we will use a linear hierarchical model
%Each trial is drawn

nTrialsPerSession = 180;
nSessions         = 4;
nParticipants     = 25;
nSimulations       = 1000; %Use at least 1000 for good estimates

%Define random distributions.
%Use function handle and param vector to make this easy to change to
%other random distribution and change params
% @normrnd - { mean standardDeviation }

%Distribution to draw each trial from
trialRandFun = @normrnd;
trialParam   = { 0 100};

%Distribution to draw each sessions bias from
sessionRandFun = @normrnd;
sessionParam = { 0 10 };

%Draw the bias from each participant from a normal distribution
%Subtle point This model means that EVERYONE has a non-zero bias, BUT the
%POPULATION MEAN is 0.  
% participantRandFun = @normrnd;
% participantParam   = { 0 1};

%Draw the bias from each participant from a normal distribution  
participantRandFun = @normrnd;
participantParam   = { 5 10};

% 
%Now let's do something tricky.  Define a function to simulate participants
%that are have a 0 bias or non-zero bias.  For those with bias draw the
%bias amount from a normal distribtion.  Heavily leveraging anonymous
%function capabilities. 
% biasMean = 0.08;
% biasStd  = 1;
% participantRandFun = @(varargin) binornd(varargin{:}) .* normrnd(biasMean,biasStd,varargin{3:end});
% participantParam   = { 1 .25}; %params are: {n p}



%Allocate a matrix to hold simulated data.
dataMatrix = NaN(nParticipants,nSessions,nTrialsPerSession);
%Holding results in a non-preallocated struct.  Can be a bit slow.
results = struct();

%
for iSim = 1:nSimulations,
    
    %clear data
    dataMatrix(:) = NaN;
    
    %Now let's simulate an experiment of data. 
    simulateSingleExperiment();
    
    %ANALYSIS BITS
    
    %First let's do the simplest analysis. This is to test the null-hypothesis
    %that the population is 0. 
    % To do thisCombine all trials/sesssions and
    
    %Use bit of a matlab trick (:,:) implicity combine last two dimensions and
    %to take a mean across both sessions and trials.
    %Could also call mean twice.
    meanData = mean(dataMatrix(:,:),2);
    % calculating mean across sessions
    meanTrialData = mean(dataMatrix,3);
    
    results(iSim).dataMatrix = dataMatrix;
    
    %Do a single-sample t-test against 0.
    [results(iSim).ttest_h results(iSim).ttest_p results(iSim).ttest_ci results(iSim).ttest_stats] = ttest(meanData);
    % Running factor analysis
    %[LAMBDA, PSI, T, STATS] = factoran(meanTrialData, 1);
    %[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(meanTrialData);
    
    [R, cp] = corrcoef(meanTrialData);
    results(iSim).R = R;
    results(iSim).cp = cp;
    
    %Cronbach's alpha - sessions
    type = 'C-k';
    [results(iSim).alpha_r, LB, UB, F, df1, df2, results(iSim).alpha_p] = ICC(meanTrialData, type);
    results(iSim).alpha_h = results(iSim).alpha_p < .05; %check if p-vale is below .05
    
    %Now let's try calculating a bayes factor.  Using a nice matlab tool
    %box.  For this we are going to compare two hypotheses.
    %Calculate the bayes factor comparing a null-hypothesis of 0 with a
    %non-zero effect using a Cauchy prior.
    results(iSim).bf = t1smpbf(results(iSim).ttest_stats.tstat,100);
    
  
end


%Now do some summary of the simulations:

%What percent of experiments found a significant effect
%This is the frequentist viewpoint.  If simulations above have a non-zero
%effect it is the power.
percentExperimentsSigTtest = mean([results(:).ttest_h])
percentExperimentsSigAlpha = mean([results(:).alpha_h])

%What is the mean bayes factor
meanBF = mean([results(:).bf])

% Average Cronbach's alpha across experiments.
meanAlpha = mean([results(:).alpha_r])

%What was the average correlation coefficient across sessions?
averageCoeff = mean([results(:).R]);



%Contain the simulation for a single experimant in a nested function. Just for making the  code
%more readable. Note: "nested" functions have full access to read and
%modify variables in the main function.  Useful, but take care not to use
%same name for different things. 
    function simulateSingleExperiment()
        
        
        for iParticipant = 1:nParticipants,
            
            %Here we use the function handle defined above to draw a bias value for
            %the responses for this participant.
            participantDraw = participantRandFun(participantParam{:});
            
            for iSession = 1:nSessions,
                
                %Here we draw a bias value to be added for the session.
                sessionDraw = sessionRandFun(sessionParam{:});
                
                %Instead of looping over each trial to draw the random
                %number just ask matlab to dtaw nTrials of random numbers
                %at once.  Turns out normrnd() has a 20microsec overhead on each
                %call that removing it from the trial loop 
                trialDraw =  trialRandFun(trialParam{:},nTrialsPerSession,1);
                dataMatrix(iParticipant,iSession,:) = trialDraw + sessionDraw + participantDraw;

                %This is just here to document a per trial loop. Removed
                %for speed. 
                %for iTrial = 1:nTrialsPerSession,                     
                %end
                
                
            end
        end
        
        
    end
end






