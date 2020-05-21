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

%% Plotting data

%% STATS - Cronbach's Alpha


%% Organising data for ANOVA