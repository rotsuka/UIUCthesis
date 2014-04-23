% inpBridge: Give the IDs of the links which you want to be bridges.
% Separate the link IDs by spaces.

function bridgeLinks=inpBridge

prompt={'Link IDs'};

% Input bridges
inp=inputdlg(prompt,'Which links are bridges?',1);

% Split link IDs
str_bridges=inp{1};
bridgeLinks=regexp(str_bridges,'\s+','split');
