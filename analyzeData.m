%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code reproduces stimuli and data analysis from the paper:
%
% Mitchell, Harris, Benstock & Ales (2020). The reliability of pseudoneglect is task dependent. Neuropsychologia
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Reset matlab state
close all
clear all




%Check if this file exists to see if we need to download data. 

dataCheckFile = fullfile('Data','P01','Analysis','Tactile','P01_TRBline1.jpg');
disp(['Checking if data has been downloaded by looking for file: ' dataCheckFile])

%Find where we are running this script from
codeLocation = mfilename('fullpath');
[codeFilePath,name,ext]=fileparts(codeLocation);

%Try to download and unpack the data if it doesn't exist. 
if ~exist(dataCheckFile,'file')
       
    url = 'https://osf.io/ch4gs/download';
    filename = 'data.zip';
    disp(['Did not find data, started downloading data from: ' url])    
    websave(fullfile(codeFilePath,filename),url);
    disp('Extracting data archive.')
    unzip(fullfile(codeFilePath,filename))
else
    disp('Found data, starting analysis');
end

%add all the code directories to the path. 
addpath(genpath(codeFilePath));

%Looking for Palamedes toolbox. 
has_pal = ~isempty(which('PAL_PFML_Fit'));
if ~has_pal
    warning('Missing Palamedes toolbox required for psychometric function fitting')
    disp('Please download the palamedes toolbox and add it to your path')
    disp('http://www.palamedestoolbox.org')
    return;
end

%Start analysis code.  

%This function changes matlab default plotting options to be a bit nicer
%looking, these changes get inherited by all subsequent plot commands. 
setPlotDefaults();

%Supress the warnings if data directories already exist. 
warning('off', 'MATLAB:MKDIR:DirectoryExists');

%This script runs the 
disp('Running the visual analysis')
pp_analysisVisual

disp('Running the tactile analysis')
pp_analysisTactile

disp('Running gorup reliability analysis')
has_finv = ~isempty(which('finv'));
if ~has_finv
    warning('Warning calculating ICC for reliability analysis requires the statitistics toolbox\n Calculating ICC will generate error')
end
%Note that there are multiple libraries that override matlab default
%statistics function and can be broken/incorrect.  If you have errors
%calculating ICC check that: finv betapdf and betainv functions are from
%matlab statistics toolbox and not a 3rd-party toolbox
%This can be done using:
%which('finv')

grp_analysisReliability


disp('Running the psychometric function fits')
pp_psychometricVisual


%Now start the analysis for the supplem ent 
disp('***********************')

disp('Now starting analysis for supplementary material')

disp('Running hand used analysis')


grp_analysis_handUsed


disp('Running line contrast analysis')
pp_psychometricContrast


grp_analysis_LMlineContrast

