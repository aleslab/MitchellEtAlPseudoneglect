# MitchellEtAlPseudoneglect
Collection of files for Investigating reliability and task demands in pseudoneglect

DOWNLOAD

Download and unzip raw Data from OSF - make sure data is saved in folder called 'Data' within repository

Download all code directly from github - place code in the same folder as 'Data' (NOT inside 'Data')


INSTRUCTIONS FOR RUNNING THE EXPERIMENT STIMULUS

To run the experiment you will need:
- MATLAB (R2017a and above) 
- The latest version of psychtoolbox (http://psychtoolbox.org/)
- To run Line bisection task you will also need a touchscreen with a Linux OS and a stylus

** Tactile rod bisection is not run through MATLAB and requires a tactile rod platform and rods **

Open 'stimulus_code' 

To run landmark task: lm.m
To run line-bisection task: mlb.m

Output files saved in same folder: 
- Landmark: LM_date_time.mat
- Line bisection: MLB_date_time.mat

All data-files are saved in individual session folders in each participant folder alongside initial analyses


INSTRUCTIONS FOR DATA ANALYSIS

To run data analysis you will need:
- MATLAB
- Palamedes toolbox to run psychophysical analyses (http://www.palamedestoolbox.org/download.html)
- Add the 'external' folder to your MATLAB path (includes necessary functions, e.g. ICC)

** NOTE: The ICC may not work on some platforms, the data to run this outside of MATLAB should be present in data **
** NOTE2: Make sure MATLAB always runs from script path - if not may get path errors **

All individual analysis files are saved so group analysis should run straight away - but individual analysis scripts are available.

~ Individual participant analyses ~

Code: folder/reliabilityanalysis_code

To analyse visual data for each participant run: pp_analysisVisual.m THEN pp_psychometricVisual.m

To analyse tactile data for each participant run: pp_analysisTactile.m

Data will end up here: folder/Data/PP/Analysis

../Visual

File: PP_visualanalysisStart.mat <- data before psychometric fitting
Results from the landmark task: lm
- average percentage of right side responses alongside raw data
Results from the line bisection task: mlb
- average mean bisection error alongside raw data

File: PP_visualanalysis <- data after psychometric fitting
No change to mlb but lm now includes PSE data for each session and overall

../Tactile

File: PP_tactileanalysis
Results: trb
- average mean bisection error alongside raw data


~ Group level analyses ~
To run analysis run: grp_analysisReliability.m
Data will end up here: folder/Analysis

grp_analysisReliability will produce 'Analysis' folder
Output: ReliabilityAnalysis.mat
- Produces a variety of different plots 
- Cronbach's alpha and t-test results
- Covariance matrices (not presented in pub)

Key structure: 'results' used to calculate alpha and make plots

Organisation of structures (by column):

Results.plotting.modalities: 
(1) observers 
(2) bisection error total mean
(3) bisection error landmark 
(4) bisection error manual line bisection
(5) bisection error tactile rod bisection
(6) bisection error standard deviation all tasks
(7) SEM all tasks
(8) observers organised by bias

Results.plotting.sessions.(task):  <- one for each lm/mlb/trb
(1) observers 
(2) bisection error task mean
(3) bisection error session 1
(4) bisection error session 2
(5) bisection error session 3
(6) bisection error session 4
(7) standard deviation task
(8) SEM task
(9) observers organised by bias


ALPHA RESULTS:

Reliability across session:
- Landmark: results.sessions.lm
- Line bisection: results.sessions.mlb
- Tactile rod: results.sessions.trb

Reliability across modalitiy:
- results.modalities.all

T-tests:
- Landmark: results.modalities.lmT
- Line bisection: results.modalities.mlbT
- Tactile rod: results.modalities.trbT


To produce figures found in text run: publicationPlots.m
Plots will end up here: folder/Analysis/Plots

This will produce 'Plots' folder within 'Analysis'
- Create and save all plots presented in publication

** 'grp_analysis_handUsed.m' will produce hand used analysis presented in supplimentary
** 'pp_psychometricContrast.m' and 'grp_analysis_LMlineContrast.m' will produce line contrast analysis presented in supplimentary
