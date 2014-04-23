% getDistanceEQ: takes the location of the earthquake and writes a new map
% object. The order here that the mapLinks are changed is not important.
%
% INPUTS
% mapLinks: a mapLinks object
% numLinks: integer for the number of links in the system
% EQ: an earthquake object

function mapLinks=getDistanceEQ(mapLinks,numLinks,EQ)

keyLinks=keys(mapLinks);

% Seperate x and y coordinates of EQ
xEQ=EQ.xycoord(1);
yEQ=EQ.xycoord(2);

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % Get x and coordinates of link centroid
    xlink=clink.xycoord(1);
    ylink=clink.xycoord(2);
    
    % Change bridges only
    if clink.isBridge==1
        
        % Distance formula
        dist=sqrt((xEQ-xlink)^2+(yEQ-ylink)^2);
        clink.distFromEQ=dist;
        
    end
    
    mapLinks(ckey)=clink;
    
end
        
        