%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Reset matlab state
close all
clear all



%Check if this file exists to see if we need to download data. 
dataCheckFile = fullfile('Data','P01','Analysis','Tactile','P01_TRBline1.jpg')

%Find where we are running this script from
codeLocation = mfilename('fullpath');
[codeFilePath,name,ext]=fileparts(codeLocation);

%Try to download and unpack the data if it doesn't exist. 
if ~exist(dataCheckFile,'file')
    
    url = 'https://osf.io/ch4gs/download'
    filename = 'data.zip'
    disp('Downloading Data')
    
    websave(fullfile(codeFilePath,filename),url);
    unzip(fullfile(codeFilePath,filename))
end

%add all the code directories to the path. 
addpath(genpath(codeFilePath));


%Start analysis code.  

%This runs the 
disp('Running the visual analysis')
pp_analysisVisual

disp('Running the psychometric function fits')
pp_psychometricVisual

disp('Running the tactile analysis')
pp_analysisTactile

disp('Running gorup reliability analysis')
grp_analysisReliability


