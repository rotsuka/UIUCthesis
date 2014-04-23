% inpLinks: Creates links by giving them an ID. State which two nodes
% connect the link (upstream downstream) and also the parameters of the
% road

function [valueLinks,keyLinks]=inpLinks(numLinks)

valueLinks=cell(1,numLinks);
keyLinks=cell(1,numLinks);

% prompt needs ID,type,conn nodes,free flow speed,backward propagating wave
% speed, maximum flow capacity, number of lanes

prompt={'ID','Type (s=start,e=end,i=internal)','Connecting node IDs',...
    'Free Flow Speed (km/hour)','Backward Propagating Wave Speed (km/hour)',...
    'Maximum Flow Cap (veh/hr/lane)','Number of Lanes'};

for i=1:numLinks
    
    inp=inputdlg(prompt,['Link ' num2str(i)],1,...
        {['l',num2str(i)],'0','node IDs','100','20','2400','2'});
    
    l=link;
    
    l.ID=inp{1};
    l.type=inp{2};
    l.nodes=inp{3};
    l.vmax=str2num(inp{4});
    l.wf=str2num(inp{5});
    l.maxCap=str2num(inp{6});
    l.noLanes=str2num(inp{7});
    
    keyLinks{i}=l.ID;
    valueLinks{i}=l;
    
end
