% getDistVec: takes in the mapLinks object and gives a column vector of the
% distances of the bridges from the EQ. The number of entries is equivalent
% to the number of bridges in the system.
%
% INPUTS
% mapLinks: a mapLinks object
% keyLinks: a cell array of link IDs
% bridgeLinks: cell array that has the IDs of links that are bridges
% numLinks: integer for the number of links in the system

function distVec=getDistVec(mapLinks,keyLinks,bridgeLinks,numLinks)

% Initialize the vector
distVec=zeros(length(bridgeLinks),1);

% Keep track of the distVec
ind=1;

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % Change bridges only
    if clink.isBridge==1
        
        % Populate and increment counter
        distVec(ind)=clink.distFromEQ;
        ind=ind+1;
        
    end

end



