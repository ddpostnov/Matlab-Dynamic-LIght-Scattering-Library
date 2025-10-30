%getEdgeSizeSLSCI - needed for fix in some of SLSCI processing steps


% Copyright @Authors: DD Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 05 May 2025

%------------- BEGIN CODE --------------

function edgeSize=getEdgeSizeSLSCI(imgK,thresh)
edgeSize=1;
while (sum(imgK(:,edgeSize)>imgK(:,edgeSize+1))/size(imgK,1))>thresh || (sum(imgK(edgeSize,:)>imgK(edgeSize+1,:))/size(imgK,2))>thresh || (sum(imgK(:,end-edgeSize+1)>imgK(:,end-edgeSize))/size(imgK,1))>thresh || (sum(imgK(end-edgeSize+1,:)>imgK(end-edgeSize,:))/size(imgK,2))>thresh
edgeSize=edgeSize+1;
end
edgeSize=edgeSize-1;
end


%------------- END OF CODE --------------