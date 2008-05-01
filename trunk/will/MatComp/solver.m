function W = solver(B)
W = zeros(0,4);

myMoves = zeros(0,5);
myConnections = zeros(0,4);
myConLen = zeros(0,1);
myID = zeros(nnz(B));

%shift everything avoid boundry conditions
thisMat = zeros(size(B)+2);
thisMat(2:end-1,2:end-1)=B;
wiredMat = true(size(B)+2);
wiredMat(2:end-1,2:end-1) = B>0;
myNumMask = zeros(size(thisMat));
myNumMask(thisMat(:)>0)=1:nnz(thisMat);

pegs = unique(nonzeros(thisMat));
numPasses=0;
maxPasses=2;

IJmat=cell2mat(arrayfun(@GetPairs,pegs,'uniformoutput',false));

while maxPasses > numPasses && ~isempty(IJmat)
    numPasses=numPasses+1;

    [bestScores inds] = sort(abs(IJmat(:,1)-IJmat(:,3)) + abs(IJmat(:,2)-IJmat(:,4)) - 2*IJmat(:,5),'ascend');

    for i = 1:nnz(bestScores<0)
        testMove = FindShortestPath(IJmat(inds(i),:));

        if ~isempty(testMove) && size(testMove,1)<2*IJmat(inds(i),5)
            wiredMat(testMove(:,1),testMove(:,2))=1;
            pos1=myNumMask(IJmat(inds(i),1),IJmat(inds(i),2));
            pos2=myNumMask(IJmat(inds(i),3),IJmat(inds(i),4));

            myConnections = [myConnections; IJmat(inds(i),1:4)];
            myConLen = [myConLen; size(testMove,1)];
            myID(pos1,pos2) = size(myConnections,1);

            myMoves=[myMoves;testMove ones(size(testMove,1),1)*size(myConnections,1)];
        end
        
        if any(sum(myID>0,2)>1)
            checks = nonzeros(myID(sum(myID>0,2)>1,:));
            %checks(checks(:)~=0)=myConLen(nonzeros(checks));
            %[vals inds] = sort(checks,2,'descend');
            %tempMask = vals>0;

            [val ind]=max(myConLen(checks));
            myConnections(checks(ind),:)=0;
            myMoves(myMoves(:,5)==checks(ind),:)=0;
            myConLen(checks(ind))=0;
            myConnections(myConnections==checks(ind))=0;
        end


    end

    
    IJmat=cell2mat(arrayfun(@GetPairs,pegs,'uniformoutput',false));
    
end

%shift everything back
W=myMoves(myMoves(:,1)~=0,1:4)-1;

    function [Moves]=FindShortestPath(Pos)
        allotedMoves=2*Pos(5);
        %since order is unimportant, always create paths with ascending
        %numbers for easy of programming

        startPos = Pos(1:2);
        endPos = Pos(3:4);

        
        Moves = zeros(allotedMoves,4);
        thiscounter = 0;
        currentPos = startPos;

        while (thiscounter<allotedMoves)
            thiscounter=thiscounter+1;
            Moves(thiscounter,1:2)=currentPos;
            currentDist=sum(abs(currentPos-endPos));
            thisMove = repmat(currentPos,[4,1])+[eye(2); -eye(2)];
            thisPoss = ~wiredMat(sub2ind(size(wiredMat),thisMove(:,1),thisMove(:,2)));
            
            thisDist = sum(abs(thisMove - repmat(endPos,[4,1])),2);

            if any(thisDist==0)
                Moves(thiscounter,3:4)=thisMove(thisDist==0,:);
                break
            end
            thisMask = (thisDist+~thisPoss.*10000)< currentDist;
            
            if ~any(thisMask) && thiscounter~=1
                thiscounter=0;
                break
            elseif ~any(thisMask) && thiscounter==1 && any(~thisPoss)

                Moves(thiscounter,3:4)=thisMove(find(~thisPoss,1),:);

            else
                Moves(thiscounter,3:4)=thisMove(find(thisMask,1),:);                

            end
            

            
            currentPos = Moves(thiscounter,3:4);
            
        end
        %if counter was reset to 0 then Moves will be empty
        Moves = Moves(1:thiscounter,:);
        return
       
    end

    function thesePairs=GetPairs(value)
        [I J]=find(thisMat==value);
        thesePairs = zeros((length(I)-1)^2,5);
        if size(I,1)<=1
            thesePairs = [];
            return
        end
        thesePairs(:,5)=value;
        counter=1;
        for k=1:length(I)
            mask=(1:length(I))>k;
            thesePairs(counter:counter+nnz(mask)-1,1:2)=repmat([I(k) J(k)],[nnz(mask) 1]);
            thesePairs(counter:counter+nnz(mask)-1,3:4)=[I(mask) J(mask)];
            counter=counter+nnz(mask)-1;
        end
        
        thesePairs=thesePairs(1:counter+1,:);
        
        alreadyParsed=ismember(thesePairs(:,1:4),myConnections,'rows');
        if any(alreadyParsed)
            thesePairs=thesePairs(~alreadyParsed,:);
        end
        
    end

end