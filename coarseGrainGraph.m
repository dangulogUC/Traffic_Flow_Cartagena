function [Q,space] = coarseGrainGraph(adj,spXY)
% This function merges intersections with degree equals to 2.


A = adj;
Aorig = A;

deg = sum((A>0),2);
aux = find(deg==2);
while(~isempty(aux))              
    for kk = 1:length(A)
        deg = sum(A(kk,:)>0);               
        if(deg==2)
            ii = find(A(kk,:)>0);
            A(ii(1),ii(2)) = Aorig(ii(1),kk)+Aorig(ii(2),kk);        
            A(ii(2),ii(1)) = Aorig(ii(1),kk)+Aorig(ii(2),kk);                
            A(kk,ii(1))=0;    
            A(ii(1),kk)=0;            
            A(kk,ii(2))=0;    
            A(ii(2),kk)=0;                
        end
    end
    deg = sum((A>0),2);
    aux = find(deg==2);
end

deg = sum((A>0),2);
nBlank1 = sum(deg==0);

G0 = graph(A);
bins = conncomp(G0);
nBlank2 = length(find(bins~=1));

nFinal = length(A)-(nBlank1+nBlank2);
space = zeros(nFinal,2);

u = 1;
for i = 1:length(A)
    if(sum(A(:,i)~=0) && bins(i)==1)
       space(u,:)=spXY(i,:);
        u = u+1;
    end    
end

degree = sum(A>0,2);
nodesElim1 = find(degree==0);
nodesElim2 = find(bins~=1);

row = [nodesElim1;nodesElim2'];        % Could be an array of rows to exclude
column = [nodesElim1;nodesElim2'];     % Could be an array of columns to exclude

Q = A(~ismember(1:size(A, 1), row), ~ismember(1:size(A, 2), column));



