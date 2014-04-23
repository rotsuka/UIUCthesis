% setInitialParams: set the jam desnity and max flow to the mapLink objects
% before the filter runs
%
% INPUTS:
% mapLinks: a mapLinks object
% numLinks: integer for the number of links in the system

function mapLinks=setInitialParams(mapLinks,numLinks)

keyLinks=keys(mapLinks);

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % Get necessary params
    wf=clink.wf;
    vmax=clink.vmax;
    
    % Get the current capacity
    clink.qmax=clink.maxCap;
    qmax=clink.qmax;
    
    % Get rhoj
    clink.rhoj=qmax/wf+qmax/vmax;
    
    % Overwrite link object
    mapLinks(ckey)=clink;
    
end
    
    