% getEQParams: gives back the magnitude and x and y locations of an EQ
% object
%
% INPUTS:
% isEQInp: integer that tells the model whether to consider the EQ or not
% trueEQ: the EQ object of the determinisitic solution
%
% NOTE: the uncertainties are hardcoded in this function

function [m x y]=getEQParams(isEQInp,trueEQ)

% Once the earthquake happens, its characteristics are known
% with a certain degree of certainty so the variances are fairly low. Thus,
% the distributions on these parameters are narrow.
sigm=.3;
sigx=5;
sigy=5;

% If there is an EQ input to the evolution model
if isEQInp==1
    
    % Draw earthquake magnitude and location from normal distributions
    % centered around the deterministic EQ object of the true solution
    m=normrnd(trueEQ.mag,sigm);
    x=normrnd(trueEQ.xycoord(1),sigx);
    y=normrnd(trueEQ.xycoord(2),sigy);

% No EQ input    
else
    
    % Force it to be insignifcant by giving such parameters that the
    % computed PGA is nearly zero
    m=0;
    x=1e6;
    y=1e6;
    
end

    