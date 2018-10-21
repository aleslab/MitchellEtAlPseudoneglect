function [results, resultsAllCond]=simulateBias(varargin)
%To simulate the data we will use a linear hierarchical model



%Define the different options we can set and their default values
p = inputParser;
addParameter(p,'nTrialsPerSession',180);
addParameter(p,'nSessions',4);
addParameter(p,'nParticipants',30);
addParameter(p,'nConditions',2);
addParameter(p,'nSimulations',100);%Use at least 1000 for good estimate
addParameter(p,'trialRandFun', @normrnd);
addParameter(p,'trialParam', {0 1; 0 1});
addParameter(p,'conditionRandFun', @normrnd);
addParameter(p,'conditionParam', {0 4});
addParameter(p,'sessionRandFun', @normrnd);
addParameter(p,'sessionParam', {0 2});
addParameter(p,'participantRandFun', @normrnd);
addParameter(p,'participantParam', {5 10});


parse(p,varargin{:});

nTrialsPerSession = p.Results.nTrialsPerSession;
nSessions         = p.Results.nSessions;
nParticipants     = p.Results.nParticipants;
nConditions       = p.Results.nConditions;
nSimulations      = p.Results.nSimulations; 

%Define random distributions.
%Use function handle and param vector to make this easy to change to
%other random distribution and change params
% @normrnd - { mean standardDeviation }

%Distribution to draw each trial from, 1 row per condition.
%This allows us to have different conditions that have different
%measurement precisions. 
trialRandFun = p.Results.trialRandFun;
trialParam   = p.Results.trialParam;

%Cross-condition variability: This allows for conditions to have different
%bias values. 
conditionRandFun = p.Results.conditionRandFun;
conditionParam =  p.Results.conditionParam;

%Distribution to draw each sessions bias from
sessionRandFun = p.Results.sessionRandFun;
sessionParam = p.Results.sessionParam;

%Draw the bias from each participant from a normal distribution
%Subtle point This model means that EVERYONE has a non-zero bias, BUT the
%POPULATION MEAN is 0.  
% participantRandFun = @normrnd;
% participantParam   = { 0 1};
  
participantRandFun = p.Results.participantRandFun;
participantParam   = p.Results.participantParam;

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
dataMatrix = NaN(nParticipants,nConditions,nSessions,nTrialsPerSession);
%Holding results in a non-preallocated struct.  Can be a bit slow.
results = struct();
resultsAllCond = struct();
%
for iSim = 1:nSimulations,
    
    %clear data
    dataMatrix(:) = NaN;
    
    %Now let's simulate an experiment of data. 
    simulateSingleExperiment();
    
    %ANALYSIS BITS
    
    %First let's do the simplest analysis. This is to test the null-hypothesis
    %that the population is 0. 
    % To do thisCombine all trials/sesssions and do the same thing for all
    % conditions

    for iCond = 1:nConditions,
        
        %Use bit of a matlab trick (:,:) implicity combine last two dimensions and
        %to take a mean across both sessions and trials.
        %Could also call mean twice.
        meanData = squeeze(mean(dataMatrix(:,iCond,:),3));
        %Calculating mean across trials for this conditions;
        meanTrialData = squeeze(mean(dataMatrix(:,iCond,:,:),4));
        
        results(iCond,iSim).dataMatrix = dataMatrix;
        
        %Do a single-sample t-test against 0.
        [results(iCond,iSim).ttest_h results(iCond,iSim).ttest_p results(iCond,iSim).ttest_ci results(iCond,iSim).ttest_stats] = ttest(meanData);
        % Running factor analysis
        %[LAMBDA, PSI, T, STATS] = factoran(meanTrialData, 1);
        %[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(meanTrialData);
        
        
        %Cronbach's alpha - sessions
        type = 'C-k';
        [results(iCond,iSim).alpha_r, LB, UB, F, df1, df2, results(iCond,iSim).alpha_p] = ICC(meanTrialData, type);
        results(iCond,iSim).alpha_h = results(iCond,iSim).alpha_p < .05; %check if p-vale is below .05
        
        %Now let's try calculating a bayes factor.  Using a nice matlab tool
        %box.  For this we are going to compare two hypotheses.
        %Calculate the bayes factor comparing a null-hypothesis of 0 with a
        %non-zero effect using a Cauchy prior.
        results(iCond,iSim).bf = t1smpbf(results(iCond,iSim).ttest_stats.tstat,100);
    end
    
    %Calculating mean across trials for this conditions;
    meanTrialData = squeeze(mean(dataMatrix(:,:,:,:),4));
    
    
    
    %Cronbach's alpha - Across all conditions/sessions. 
        type = 'C-k';
        [resultsAllCond(iSim).all_alpha_r, ~, ~, ~, ~, ~, resultsAllCond(iSim).all_alpha_p] = ICC(meanTrialData(:,:), type);
        resultsAllCond(iSim).all_alpha_h = resultsAllCond(iSim).all_alpha_p < .05; %check if p-vale is below .05
        
     %Next let's take average across sessions and correlate   
     meanSessionData = squeeze(mean(meanTrialData,3));
         [resultsAllCond(iSim).session_mean_alpha_r, ~, ~, ~, ~, ~, resultsAllCond(iSim).session_mean_alpha_p] = ICC(meanSessionData(:,:), type);
          resultsAllCond(iSim).session_mean_alpha_h = resultsAllCond(iSim).session_mean_alpha_p < .05; %check if p-vale is below .05

                    
     %Correlation Coefficient across all sessions. 
         [R, cp] = corrcoef(meanTrialData(:,:));
         resultsAllCond(iSim).all_corrcoeff = R;
         resultsAllCond(iSim).all_corrcoeff_p = cp;

          
     %Correlation Coefficient across conditions. 
         [R, cp] = corrcoef(meanSessionData);
         resultsAllCond(iSim).session_mean_corrcoeff = R;
         resultsAllCond(iSim).session_mean_corrcoeff_p = cp;

  
end


for iCond = 1:nConditions,
%Now do some summary of the simulations:

%What percent of experiments found a significant effect
%This is the frequentist viewpoint.  If simulations above have a non-zero
%effect it is the power.
percentExperimentsSigTtest(iCond) = mean([results(iCond,:).ttest_h],2);
percentExperimentsSigAlpha(iCond) = mean([results(iCond,:).alpha_h],2);

%What is the mean bayes factor
meanBF(iCond) = mean([results(iCond,:).bf],2);

% Average Cronbach's alpha across sessions?.
meanAlpha(iCond) = mean([results(iCond,:).alpha_r],2);

%What was the average correlation coefficient across sessions?
%averageCoeff(iCond) = mean([results(iCond,:).R],2);

disp(['Results for condition: ' num2str(iCond)]);
disp(['Proporation of inferential tests significant: ']);
disp(['Ttest: ' num2str(percentExperimentsSigTtest(iCond)) ...
    ' Cronbach''s Alpha: ' num2str(percentExperimentsSigAlpha(iCond))]);
disp(['Mean Cronbach''s Alpha: ' num2str(meanAlpha(iCond))]);
disp('-------------------')

end

disp('Results combining conditions, keeping sessions separate')
disp(['Proportion significant: ' num2str(mean([resultsAllCond(:).all_alpha_h]))])
disp(['Mean Cronbach''s Alpha: ' num2str(mean([resultsAllCond(:).all_alpha_r]))]);

disp('Results combining conditions, and averaging sessions together')
disp(['Proportion significant: ' num2str(mean([resultsAllCond(:).session_mean_alpha_h]))])
disp(['Mean Cronbach''s Alpha: ' num2str(mean([resultsAllCond(:).session_mean_alpha_r]))]);


meanCondCorrCoeff= squeeze(mean(cat(3,resultsAllCond(:).session_mean_corrcoeff),3));
meanSessCorrCoeff= squeeze(mean(cat(3,resultsAllCond(:).all_corrcoeff),3));
disp('----------')
disp('Correlation Coefficients across conditions')
disp(meanCondCorrCoeff)
disp('----------')
disp(meanSessCorrCoeff)
disp('----------')

%Contain the simulation for a single experimant in a nested function. Just for making the  code
%more readable. Note: "nested" functions have full access to read and
%modify variables in the main function.  Useful, but take care not to use
%same name for different things. 
    function simulateSingleExperiment()
        
        
        for iParticipant = 1:nParticipants,
            
            %Here we use the function handle defined above to draw a bias value for
            %the responses for this participant.
            participantDraw = participantRandFun(participantParam{:});
            
            for iCondition = 1:nConditions,
                
                %Here we draw a a cross-condition variability
                conditionDraw   = conditionRandFun(conditionParam{:});
                
                for iSession = 1:nSessions,
                    
                    %Here we draw a bias value to be added for the session.
                    sessionDraw = sessionRandFun(sessionParam{:});
                    
                    %Instead of looping over each trial to draw the random
                    %number just ask matlab to dtaw nTrials of random numbers
                    %at once.  Turns out normrnd() has a 20microsec overhead on each
                    %call that removing it from the trial loop
                    trialDraw =  trialRandFun(trialParam{iCondition,:},nTrialsPerSession,1);
                    dataMatrix(iParticipant,iCondition,iSession,:) = trialDraw + sessionDraw + participantDraw + conditionDraw;
                    
                    %This is just here to document a per trial loop. Removed
                    %for speed.
                    %for iTrial = 1:nTrialsPerSession,
                    %end
                    
                end
                
                
            end
        end
        
        
    end
end






