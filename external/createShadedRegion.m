function [ hLine, hPatch ] = createShadedRegion( x,y,yLo,yHi,varargin )
%createShadedPlot Plots a line with a shaded region around it. 
%   
%[ hLine hPatch ] = createShadedPlot( x,y,yLo,yHi, [extra inputs passed to plot] )
%
%This is an almost drop in replacement for plot that plots with a shaded
%region around the line.  
%
%4 inputs required:
%
% x, x location to plot. 
% y- y values to plot and location x
% yLo - Lower value of the shaded region to plot
% yHi - Upper location of the shaded region. 
%
%Any extra inputs are passed to the plot() function for formating the line:
%
%Example:
%x = linspace(0,1,100);
%y = x*2;
%createShadedRegion(x,y,y+2,y-2,':','color',[1 0 0])

%TODO: Error checking!

%These points will be used to create the patch that we will make
%transparent
px = [x(:); flipud(x(:))];    % use px=[x1;x2(end:-1:1)]; if row vecs
py = [ yLo(:) ; flipud(yHi(:)) ];    % use ...

%Plot the line. 
hLine = plot(x,y,varargin{:});

hPatch=patch(px,py,hLine.Color,'linestyle','none','facealpha',.25);

end
