%% AG.Mitchell 30.04.20
% Script to investigate whether hand used to complete MLB and TRB task
% affects bisection error
% And whether reliability across time is increased if error is from one
% hand only 
% Addressing reviewer 2 comments in Neuropsychologia re-submission (1.1)

%% Setting paths
clear all
nParticipants = [1:19,21:24,26:30];
% Directory
dirBias = ('M:\Alex_Files\Experiments\Bias'); %subject to change depending on where you analyse
dirAna = ('M:\Alex_Files\Experiments\Bias\Analysis');

for p = 1:length(nParticipants)
    ppID = sprintf('P%0*d',2,nParticipants(p)); %for use when navigating files
    visualFilename = sprintf('%s_visualanalysisStart.mat', ppID);
    tactileFilename = sprintf('%s_tactileanalysis.mat', ppID);
    nSessions = 1:4; %vector number of sessions each participant does
    
    % individual participant directory, analysis
    dirPP = [dirBias filesep 'Data' filesep ppID filesep 'Analysis']; %participant directory
    dirVis = [dirPP filesep 'Visual'];
    dirTac = [dirPP filesep 'Tactile'];
    cd(dirPP)
    
%% Extracting data
    % for MLB - go to visual folder
    cd(dirVis)
    load(visualFilename)
    cd(dirTac)
    load(tactileFilename)
    
    %% Getting data per session, compiling into matrices
    for s = 1:length(nSessions)
        cd(dirPP)
        session = sprintf('Session%0*d',2,nSessions(s));
        % MLB
        % left hand
        res_mlb.left_hand.error(p,s) = (mlb.(sprintf('%s', session)).error.left_hand(1))*10; %mean
        res_mlb.left_hand.sd(p,s) = (mlb.(sprintf('%s', session)).error.left_hand(2))*10; %sd
        % right hand
        res_mlb.right_hand.error(p,s) = (mlb.(sprintf('%s', session)).error.right_hand(1))*10; %mean
        res_mlb.right_hand.sd(p,s) = (mlb.(sprintf('%s', session)).error.right_hand(2))*10; %sd
        
        % TRB
        % left hand
        res_trb.left_hand.error(p,s) = (trb.(sprintf('%s', session)).error.left_hand(1))*10; %mean
        res_trb.left_hand.sd(p,s) = (trb.(sprintf('%s', session)).error.left_hand(2))*10; %sd
        % right hand
        res_trb.right_hand.error(p,s) = (trb.(sprintf('%s', session)).error.right_hand(1))*10; %mean
        res_trb.right_hand.sd(p,s) = (trb.(sprintf('%s', session)).error.right_hand(2))*10; %sd       
    end

end

%% Outlier removal
% Any outliers identified through MATLAB - more than three scaled mean
% absolute deviations away
% Loading outlier from data-set so the outliers are the same across analyses
cd(dirAna)
load('outliers.mat')
remove = find(outlierSum > 0.1); %identifies paticipants that need removing
remove = sort(remove); %sorting so in participant order
% idenitfying outliers, making NaN
for r = 1:length(remove)
    % MLB
    res_mlb.left_hand.error(remove(r)) = NaN;
    res_mlb.left_hand.sd(remove(r)) = NaN;
    res_mlb.right_hand.error(remove(r)) = NaN;
    res_mlb.right_hand.sd(remove(r)) = NaN;
    % TRB
    res_trb.left_hand.error(remove(r)) = NaN;
    res_trb.left_hand.sd(remove(r)) = NaN;
    res_trb.right_hand.error(remove(r)) = NaN;
    res_trb.right_hand.sd(remove(r)) = NaN;
end
% then removing and saving to different matrix
% MLB
res_mlb.left_hand.Rerror = res_mlb.left_hand.error(isfinite(res_mlb.left_hand.error(:,1)),:);
res_mlb.left_hand.Rsd = res_mlb.left_hand.sd(isfinite(res_mlb.left_hand.sd(:,1)),:);
res_mlb.right_hand.Rerror = res_mlb.right_hand.error(isfinite(res_mlb.right_hand.error(:,1)),:);
res_mlb.right_hand.Rsd = res_mlb.right_hand.sd(isfinite(res_mlb.right_hand.sd(:,1)),:);
% TRB
res_trb.left_hand.Rerror = res_trb.left_hand.error(isfinite(res_trb.left_hand.error(:,1)),:);
res_trb.left_hand.Rsd = res_trb.left_hand.sd(isfinite(res_trb.left_hand.sd(:,1)),:);
res_trb.right_hand.Rerror = res_trb.right_hand.error(isfinite(res_trb.right_hand.error(:,1)),:);
res_trb.right_hand.Rsd = res_trb.right_hand.sd(isfinite(res_trb.right_hand.sd(:,1)),:);

%% Averaging across session
% MLB, left hand
res_mlb.left_hand.average(:,1) = nanmean(res_mlb.left_hand.Rerror, 2); %mean
res_mlb.left_hand.average(:,2) = nanstd(res_mlb.left_hand.Rerror, 0, 2); %std
% right hand
res_mlb.right_hand.average(:,1) = nanmean(res_mlb.right_hand.Rerror, 2); %mean
res_mlb.right_hand.average(:,2) = nanstd(res_mlb.right_hand.Rerror, 0, 2); %std
% TRB
res_trb.left_hand.average(:,1) = nanmean(res_trb.left_hand.Rerror, 2); %mean
res_trb.left_hand.average(:,2) = nanstd(res_trb.left_hand.Rerror, 0, 2); %std
res_trb.right_hand.average(:,1) = nanmean(res_trb.right_hand.Rerror, 2); %mean
res_trb.right_hand.average(:,2) = nanstd(res_trb.right_hand.Rerror, 0, 2); %std

% Average across hand for each task - for plotting
mlball = [res_mlb.left_hand.average(:,1), res_mlb.right_hand.average(:,1)];
res_mlb.average = mean(mlball, 2);
trball = [res_trb.left_hand.average(:,1), res_trb.right_hand.average(:,1)];
res_trb.average = mean(trball, 2);

% Compile averages into modalities structure to use for analysis
% Load and add the landmark data here too - for comparison
load('ReliabilityAnalysis.mat', 'results')
landmark = results.plotting.modalities(:,3);
res_mods.left_hand = [res_mlb.left_hand.average(:,1), res_trb.left_hand.average(:,1), landmark];
res_mods.right_hand = [res_mlb.right_hand.average(:,1), res_trb.right_hand.average(:,1), landmark];
% data for each hand, in same modality
res_mlb.hand = [res_mlb.right_hand.average(:,1), res_mlb.left_hand.average(:,1)];
res_trb.hand = [res_trb.right_hand.average(:,1), res_trb.left_hand.average(:,1)];

%% STATS - Cronbach's Alpha
% Using Cronbach's alpha to assess reliability of mlb and trb across
% sessions for each hand (left and right) to see if splitting across hand
% used improves reliability of results (reviewer request)
type = 'C-k'; %type of ICC used - 2-k, 2-way fixed effects
% Reliability across session, for each hand
% MLB, left hand
[res_mlb.reliability.left_hand.r, res_mlb.reliability.left_hand.bound(1), res_mlb.reliability.left_hand.bound(2), ...
    res_mlb.reliability.left_hand.F, res_mlb.reliability.left_hand.df(1), res_mlb.reliability.left_hand.df(2), ...
    res_mlb.reliability.left_hand.p] = ICC(res_mlb.left_hand.Rerror, type);
% MLB, right hand
[res_mlb.reliability.right_hand.r, res_mlb.reliability.right_hand.bound(1), res_mlb.reliability.right_hand.bound(2), ...
    res_mlb.reliability.right_hand.F, res_mlb.reliability.right_hand.df(1), res_mlb.reliability.right_hand.df(2), ...
    res_mlb.reliability.right_hand.p] = ICC(res_mlb.right_hand.Rerror, type);
% TRB, left hand
[res_trb.reliability.left_hand.r, res_trb.reliability.left_hand.bound(1), res_trb.reliability.left_hand.bound(2), ...
    res_trb.reliability.left_hand.F, res_trb.reliability.left_hand.df(1), res_trb.reliability.left_hand.df(2), ...
    res_trb.reliability.left_hand.p] = ICC(res_trb.left_hand.Rerror, type);
% TRB, right hand
[res_trb.reliability.right_hand.r, res_trb.reliability.right_hand.bound(1), res_trb.reliability.right_hand.bound(2), ...
    res_trb.reliability.right_hand.F, res_trb.reliability.right_hand.df(1), res_trb.reliability.right_hand.df(2), ...
    res_trb.reliability.right_hand.p] = ICC(res_trb.right_hand.Rerror, type);
    
% Reliability across modality, for each hand - including landmarks for both (hand
% used did not differentiate)
% left hand
[res_mods.reliability.left_hand.r, res_mods.reliability.left_hand.bound(1), res_mods.reliability.left_hand.bound(2), ...
    res_mods.reliability.left_hand.F, res_mods.reliability.left_hand.df(1), res_mods.reliability.left_hand.df(2), ...
    res_mods.reliability.left_hand.p] = ICC(res_mods.left_hand, type);
% right hand
[res_mods.reliability.right_hand.r, res_mods.reliability.right_hand.bound(1), res_mods.reliability.right_hand.bound(2), ...
    res_mods.reliability.right_hand.F, res_mods.reliability.right_hand.df(1), res_mods.reliability.right_hand.df(2), ...
    res_mods.reliability.right_hand.p] = ICC(res_mods.right_hand, type);

% Correlate left and right hand values between mlb and trb
[rho,pval] = corr(res_mods.left_hand(:,1), res_mods.left_hand(:,2));
res_mods.corr.left_hand.r = rho;
res_mods.corr.left_hand.p = pval;
[rho,pval] = corr(res_mods.right_hand(:,1), res_mods.right_hand(:,2));
res_mods.corr.right_hand.r = rho;
res_mods.corr.right_hand.p = pval;

%% Cronbach's alpha across hands
% MLB
[res_mlb.reliability.hand.r, res_mlb.reliability.hand.bound(1), res_mlb.reliability.hand.bound(2), ...
    res_mlb.reliability.hand.F, res_mlb.reliability.left_hand.df(1), res_mlb.reliability.hand.df(2), ...
    res_mlb.reliability.hand.p] = ICC(res_mlb.hand, type);
% TRB
[res_trb.reliability.hand.r, res_trb.reliability.hand.bound(1), res_trb.reliability.hand.bound(2), ...
    res_trb.reliability.hand.F, res_trb.reliability.left_hand.df(1), res_trb.reliability.hand.df(2), ...
    res_trb.reliability.hand.p] = ICC(res_trb.hand, type);

%% Plotting data
plotting = struct;
plotting.observers = 1:length(res_mlb.average);

%% Organising data for ANOVA

%% saving
cd(dirAna)
save('hand_analysis.mat', 'res_mlb', 'res_trb', 'res_mods');