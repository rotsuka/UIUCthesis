classdef link
    % link: represents a stretch of road in the CTM network
    
    properties
        
        % Local to link
        ID
        type
        nodes
        xycoord=[0 0];
        roadLength=0;
        rhoj=0; % per lane
        vmax=0;
        wf=0;
        rho0=0; % total
        maxCap=0; % per lane
        qmax=0; % per lane
        
        sensorObjs % tells what sensor IDs are on each link (added later)
     
        dx=0; % computed later in code
        noCells=0; % computed later in code

        % Globally
        startCell
        endCell
        
        % These represent the flows out of the links to connecting links
        % NOT the flows within the link itself
        trueInflow=0;
        trueOutflow=0;
        
        % Paremeters for merges and diverges. Only gets computed/used if
        % link is part of a merge/diverge
        noLanes=0;
        alpha=0;
        beta=0;
        
        % Variable that says whether or not it is actually a bridge
        isBridge=0; % 0=no,1=yes
        distFromEQ=0;
      
    end
    
    methods
    end
    
end

