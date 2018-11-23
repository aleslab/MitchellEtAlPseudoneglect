function [errorbarHandle] = plot2afc( x, nCorrect, nTrials )
%plot2afc plots 2afc data with 95% confidence intervals
%   Detailed explanation goes here

percentCorrect = nCorrect./nTrials;
lowerCi = percentCorrect - binoinv(.025,nTrials,percentCorrect)./nTrials;
upperCi = binoinv(.975,nTrials,percentCorrect)./nTrials - percentCorrect;

errorbarHandle=errorbar(x,percentCorrect,lowerCi,upperCi,'o','markersize',10);
box(gca,'off');  

end
