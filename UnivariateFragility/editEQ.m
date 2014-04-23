% editEQ: creates a new EQ object by overwriting the old one.
%
% INPUTS
% m: magnitude
% xc: location in the x coordinate
% yc: location in the y coordinate

function eqObj=editEQ(m,xc,yc)

% New object
eqObj=EQ;

% Assign characteristics
eqObj.mag=m;
eqObj.xycoord=[xc yc];
