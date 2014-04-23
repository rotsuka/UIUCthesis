% inpEQ: Create the earthquake object by specifying the magnitude and
% [x y] coordinates.

function newEQ=inpEQ

prompt={'Magnitude','XY Coordinates'};

% Input bridges
inp=inputdlg(prompt,'Which links are bridges?',1,{'0','0 0'});

newEQ=EQ;

newEQ.mag=str2num(inp{1});
newEQ.xycoord=str2num(inp{2});