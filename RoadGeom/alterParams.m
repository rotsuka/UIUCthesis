% alterParams

function nvalueLinks=alterParams(valueLinks)

% Initialization
nvalueLinks=valueLinks;
numLinks=length(valueLinks);

% Noise (std dev %)
maxCapNoise=5;
vmaxNoise=5;
wfNoise=5;

for i=1:numLinks
    
    % Get current link
    nlink=valueLinks{i};
    
    % Alter parameters
    nlink.maxCap=normrnd(nlink.maxCap,(maxCapNoise/100)*nlink.maxCap);
    nlink.vmax=normrnd(nlink.vmax,(vmaxNoise/100)*nlink.vmax);
    nlink.wf=normrnd(nlink.wf,(wfNoise/100)*nlink.wf);
    
    % Replace link
    nvalueLinks{i}=nlink;
    
end

        
        

