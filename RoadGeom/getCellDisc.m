% getCellDisc: Compute dx,number of cells, and start/end cell for link
% using CFL condition for the true model
%
% INPUTS
% delt: Delta T (time step)
% mapLinks: a mapLinks object
% keyLinks: a cell array of link IDs
% numLinks: integer for the number of links

function mapLinks=getCellDisc(delt,mapLinks,keyLinks,numLinks)

% Keep track of the total number of cells
count=0;

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % Internal links
    if strcmp(clink.type,'i')==1
    
        % Compute dx based on the CFL condition: 1>=(DT/dx)*v_max
        clink.dx=(delt/3600)*clink.vmax; % discrete space step [km]
        
        % Compute the number of cells for each link
        clink.noCells=floor(clink.roadLength/...
            clink.dx); % cells, rounded down
        
        % Recalculate dx based on the rounded down cell number. This deltax
        % should be higher than the previous one, thus still satisfying the CFL
        % condition
        clink.dx=clink.roadLength/clink.noCells;
        
        % Determine start and end cell
        clink.startCell=count+1;
        clink.endCell=count+clink.noCells;
        
    else % boundary link
        
        clink.noCells=1;
        
        clink.startCell=count+1;
        clink.endCell=count+1;
        
    end
    
    % Update the count
    count=count+clink.noCells;
    
    mapLinks(ckey)=clink;
    
end