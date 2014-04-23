% getMapParamsN: overwrites the old mapLinks object with the damaged
% mapLinks object. THIS IS A SPECIAL INSTANCE OF THE FUNCTION WHERE WE
% WANT THE BRIDGE TO BE DAMAGED!
%
% INPUTS
% mapLinks: a mapLinks object
% numLinks: integer for the number of links in the system
% keyLinks: a cell array of link IDs
% DSlims: vector designating ranges of different damage states
%
% NOTE: In the true model, we force a damage state, depend on the type of
% earthquake object created (i.e. it is deterministic)

function mapLinks=getMapParamsN(mapLinks,numLinks,keyLinks,DSlims)

% Keep track of the number of bridges
ind=1;

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % Get the max capacity
    maxCap=clink.maxCap;
    
    % Check to see if the current link is a bridge
    if clink.isBridge==1
        
        isDamaged=0;
        
        while isDamaged==0
            
            % Assign preferred damage to clink 
            clink.qmax=0*maxCap;
            
            % Display text
            disp(['Bridge ' num2str(ind) ', Capacity: ' num2str(clink.qmax)]);
            
            if clink.qmax<maxCap
                
                isDamaged=1;
                
            end
            
        end

        % clink.qmax=maxCap;
        
        % Increment the index
        ind=ind+1;
        
    else % not a bridge, and thus no reduction in capacity
        
        clink.qmax=maxCap;
        
    end
    
    % Overwrite old link object
    mapLinks(ckey)=clink;
    
end
