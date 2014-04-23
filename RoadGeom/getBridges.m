% getBridges: determines the links that are bridges and overwrites map
% object
%
% INPUTS
% mapLinks: a mapLinks object
% numLinks: integer for the number of links
% bridgeLinks: cell array that has the IDs of links that are bridges

function mapLinks=getBridges(mapLinks,numLinks,bridgeLinks)

% Get keys
keyLinks=keys(mapLinks);

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    for j=1:length(bridgeLinks)
        
        if strcmp(clink.ID,bridgeLinks{j})==1
            
            clink.isBridge=1;
            
        end
        
    end
    
    % Rewrite link object
    mapLinks(ckey)=clink;
    
end