function W = solver(B)
W = zeros(0,4);

%shift everything avoid boundry conditions
thisMat = zeros(size(B)+2);
thisMat(2:end-1,2:end-1)=B;
wiredMat = true(size(B)+2);
wiredMat(2:end-1,2:end-1) = B>0;



pegs = unique(nonzeros(thisMat));

IJpairs = cell(length(pegs),1);
for i = 1:length(pegs)
    [I J]=find(thisMat==pegs(i));
    thisGroup = zeros((length(I)-1)^2,5);
    if size(thisGroup,1)==0
        continue
    end
    thisGroup(:,5)=pegs(i);
    counter=1;
    for k=1:length(I)
        mask=(1:length(I))>k;
        thisGroup(counter:counter+nnz(mask)-1,1:2)=repmat([I(k) J(k)],[nnz(mask) 1]);
        thisGroup(counter:counter+nnz(mask)-1,3:4)=[I(mask) J(mask)];
        counter=counter+nnz(mask)-1;
    end
    IJpairs{i} = thisGroup(1:counter+1,:);
end

IJmat = cell2mat(IJpairs);

[bestScores inds] = sort(abs(IJmat(:,1)-IJmat(:,3)) + abs(IJmat(:,2)-IJmat(:,4)) - 2*IJmat(:,5),'ascend');


for i = 1:nnz(bestScores<0)
    thisMove = FindShortestPath(IJmat(inds(i),:));
    
    if ~isempty(thisMove) && size(thisMove,1)<2*IJmat(inds(i),5)
        wiredMat(thisMove(:,1),thisMove(:,2))=1;
        W=[W;thisMove];
    end
end

%shift everything back
W=W-1;

    function [Moves]=FindShortestPath(Pos)
        allotedMoves=2*Pos(5);
        %since order is unimportant, always create paths with ascending
        %numbers for easy of programming
        if Pos(1) < Pos(3) || Pos(1) == Pos(3) && Pos(2) < Pos(4)
            startPos = Pos(1:2);
            endPos = Pos(3:4);
        else
            startPos = Pos(3:4);
            endPos = Pos(1:2);
        end
        
        Moves = zeros(allotedMoves,4);
        thiscounter = 0;
        currentPos = startPos;

        while (counter<allotedMoves)
            thiscounter=thiscounter+1;
            Moves(thiscounter,1:2)=currentPos;
            currentDist=sum(abs(currentPos-endPos));
            thisMove = repmat(currentPos,[2,1])+eye(2);
            thisDist = sum(abs(thisMove - repmat(endPos,[2,1])),2);
            if any(thisDist==0)
                Moves(thiscounter,3:4)=thisMove(thisDist==0,:);
                break
            end
            if (~wiredMat(thisMove(1,1),thisMove(1,2)) && thisMat(thisMove(1,1),thisMove(1,2)) ~= Pos(5)) && thisDist(1)<currentDist
                Moves(thiscounter,3:4)=thisMove(1,:);
            elseif (~wiredMat(thisMove(2,1),thisMove(2,2)) && thisMat(thisMove(2,1),thisMove(2,2)) ~= Pos(5)) && thisDist(2)<currentDist
                Moves(thiscounter,3:4)=thisMove(2,:);
            else
                thiscounter=0;
                break
            end
            
            currentPos = Moves(thiscounter,3:4);
            
        end
        %if counter was reset to 0 then Moves will be empty
        Moves = Moves(1:thiscounter,:);
        return
        
        
    end

end