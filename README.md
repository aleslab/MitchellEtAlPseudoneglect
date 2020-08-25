# MitchellEtAlPseudoneglect
Collection of files for Investigating reliability and task demands in pseudoneglect

DOWNLOAD
Download and unzip raw Data from OSF - make sure data is saved in folder called 'Data' within repository
Download all code directly from github - place code in the same folder as 'Data' (NOT inside 'Data')


INSTRUCTIONS FOR RUNNING THE EXPERIMENT STIMULUS
To run the experiment you will need:
- MATLAB (R2017a and above) 
- The latest version of psychtoolbox (http://psychtoolbox.org/)


INSTRUCTIONS FOR DATA ANALYSIS
To run data analysis you will need:
- MATLAB
- Palamedes toolbox to run psychophysical analyses (http://www.palamedestoolbox.org/download.html)
- Add the 'external' folder to your MATLAB path (includes necessary functions, e.g. ICC)

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
To produce figures found in text run: publicationPlots.m

This code will produce figures