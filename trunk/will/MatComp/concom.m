function [score W Wlab C B] = concom(B,W)
% CONCOM Finds connected components and scores the solution
%      Inputs:
%         B - Board, a matrix with component ID's, 0's means no component.
%         W - Wiring, a matrix (nx4) matrix with wires, invalid wires are
%             ignored but charged into the score. 
%
% The MATLAB Contest Team 
% April, 2008 

[ro,co] = size(B);
dW = abs(W(:,[1 2])-W(:,[3 4]));
% Extract bridges from W
Bridges = W(~any(dW,2),[1 2]);
% Calculate the wiring cost: (Bridge cost is 25, but 1 is alredy counted in W
wiringCost = size(W,1) + 24 * size(Bridges,1); 
% Keep only valid wires in W
W = W(all(W(:,[1 3])>0 & W(:,[1 3])<=ro & W(:,[2 4])>0 & W(:,[2 4])<=co,2) & (sum(dW,2)==1),:);
% Change wiring coordinates to linear indexing:
V = sort([sub2ind([ro,co],W(:,1),W(:,2)) sub2ind([ro,co],W(:,3),W(:,4))],2);
BridgeIdx = sub2ind([ro,co],Bridges(:,1),Bridges(:,2));
% Create a graph (E (edges)):
E = zeros(size(V));
[node2idx,dummy,E(:)]= unique(V);
nN = numel(node2idx);
idx2node(node2idx) = 1:nN;
% Build an adjacency list for edges:
A = zeros(nN,4);  % [v ^ > <]
for i = 1:nN
    [j,k] = find(node2idx(i)==V);
    A(i,2*(diff(V(j,:),[],2)>1)+k) = sum(E(j,:),2)-i;
end
% Keep bridges that exist over wires and are not over components
BridgeIdx = intersect(BridgeIdx,V(:));
BridgeIdx = BridgeIdx(B(BridgeIdx)==0);
% Remove bridges that do not separate vertical and horizontal wires
BridgeIdx = BridgeIdx(any(A(idx2node(BridgeIdx),[1 2]),2) & any(A(idx2node(BridgeIdx),[3 4]),2));
% Add extra nodes where Bridges exist:
A = [A;zeros(numel(BridgeIdx),4)];
node2idx = [node2idx; BridgeIdx];
j = nN;
for i = idx2node(BridgeIdx)
    j = j+1;
    A(j,:) = A(i,:);
    if A(i,3) 
        A(A(i,3),4)=j; 
    end
    if A(i,4) 
        A(A(i,4),3)=j; 
    end
    A(i,[3 4]) = 0;
    A(j,[1 2]) = 0;
end
% Label nodes
labels = zeros(size(A,1),1);
tv = find(~labels,1);
cc = 1;
while ~isempty(tv)
    while ~isempty(tv)
        labels(tv) = cc;
        tv = nonzeros(A(tv,:));
        tv = tv(labels(tv) == 0);
    end
    cc = cc+1;
    tv = find(~labels,1);
end

% Copy labels to a board
C = zeros(ro,co);
C(node2idx) = labels;
H = ~C & B>0;
C(H) = cc+(0:nnz(H)-1);    
C(BridgeIdx)=-1;

% Label wires
Wlab = zeros(size(V,1),1);
for i=1:size(A,1)
    for j = 1:4
        if A(i,j)
            Wlab(V(:,1)==node2idx(i) & V(:,2)==node2idx(A(i,j))) = labels(i);
        end
    end
end

% find short-circuited components and remove them:
P = unique(nonzeros(B));
for i = 1:cc-1
    p = unique(nonzeros(B(C==i))); 
    if numel(p)>1
        P = setdiff(P,p);
    end
end
% consider the largest island for every non-short-circuited component and
% remove it from the board:
for i = 1:numel(P)
    [sz,ch] = max(sparse(C(B==P(i)),1,1));
    if sz>1
       B(C==ch)=0;
    end
end

score = sum(nonzeros(B)) + wiringCost;
    
