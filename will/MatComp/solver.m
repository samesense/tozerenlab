function W = solver(B)
W = zeros(0,4);

myMoves = zeros(0,5);
myConnections = zeros(0,7);
myConLen = zeros(0,1);
myID = zeros(nnz(B));

MovingMat=[eye(2); -eye(2)];

%shift everything avoid boundry conditions
thisMat = zeros(size(B)+2);
thisMat(2:end-1,2:end-1)=B;
wiredMat = true(size(B)+2);
wiredMat(2:end-1,2:end-1) = B>0;
labeledWire = double(wiredMat);
myNumMask = zeros(size(thisMat));
myNumMask(thisMat(:)>0)=1:nnz(thisMat);

pegs = unique(nonzeros(thisMat));
numPasses=0;
maxPasses=4;

IJmat=cell2mat(arrayfun(@GetPairs,pegs,'uniformoutput',false));
% extragets = 0;

while ~isempty(IJmat)
    numPasses=numPasses+1;

    [bestScores inds] = sort(abs(IJmat(:,1)-IJmat(:,3)) + abs(IJmat(:,2)-IJmat(:,4)) - 2*IJmat(:,5),'ascend');

    for i = 1:nnz(bestScores<0)
        testMove = FindShortestPath(IJmat(inds(i),:),1);

        if ~isempty(testMove) && size(testMove,1)<2*IJmat(inds(i),5)
            
            wiredMat((testMove(:,2)-1)*size(wiredMat,1)+testMove(:,1))=1;
            labeledWire((testMove(:,2)-1)*size(wiredMat,1)+testMove(:,1))=IJmat(inds(i),5);
            pos1=myNumMask(IJmat(inds(i),1),IJmat(inds(i),2));
            pos2=myNumMask(IJmat(inds(i),3),IJmat(inds(i),4));

            if sum(myID(pos1,:)~=0)<1
                myConnections = [myConnections; IJmat(inds(i),1:4) IJmat(inds(i),5) pos1 pos2];
                myConLen = [myConLen; size(testMove,1)];
                myID(pos1,pos2) = size(myConnections,1);

                myMoves=[myMoves;testMove ones(size(testMove,1),1)*size(myConnections,1)];

                
            end
            

            

        end

    end

%     for k=1:nnz(B)
%         if nnz(myConnections(:,6:7)==k) > 2
%             AllConnects=find(sum(myConnections(:,6:7)==k,2)>0);
%             
%             OtherPegs=nonzeros(myConnections(AllConnects,6:7).*(myConnections(AllConnects,6:7)~=k));
%             
%             NumConnects=sum(myID(OtherPegs,:)>0,2)+sum(myID(:,OtherPegs)>0,1)';
%             
%             toLose = NumConnects<2;
%             
%             removeConnects=[nonzeros(myID(k,OtherPegs(toLose)));nonzeros(myID(OtherPegs(toLose),k))];
%             
%             %visualize(B,myMoves(:,1:4)-1);
%             %visualize(B,myMoves(ismember(myMoves(:,5),removeConnects),1:4)-1);
%             
%             myConnections(removeConnects,:)=0;
%             
%             myMoves(ismember(myMoves(:,5),removeConnects),:)=0;
%             %visualize(B,myMoves(myMoves(:,1)~=0,1:4)-1);
%         end
%         
%     end


    allCons=[myMoves(myMoves(:,1)~=0&myMoves(:,4)~=0,1:2);myMoves(myMoves(:,1)~=0&myMoves(:,4)~=0,3:4)];
    wiredMat(2:end-1,2:end-1) = B>0;

    wiredMat((allCons(:,2)-1)*size(wiredMat,1)+allCons(:,1))=1;

    if numPasses > maxPasses
        break
    end

    MovingMat=MovingMat(randperm(4),:);

    IJmat=cell2mat(arrayfun(@GetPairs,pegs,'uniformoutput',false));

end

% extragets
%shift everything back
W=unique(myMoves(myMoves(:,1)~=0&myMoves(:,4)~=0,1:4)-1,'rows');

    function [Moves]=FindShortestPath(Pos,order)
        allotedMoves=2*Pos(5);
        %since order is unimportant, always create paths with ascending
        %numbers for easy of programming
        thisWired = wiredMat - (labeledWire==Pos(end));
        wiredInds=find((labeledWire==Pos(end)));
        if order==1
            startPos = Pos(1:2);
            endPos = Pos(3:4);
        elseif order==2
            startPos = Pos(3:4);
            endPos = Pos(1:2);
        else
            Moves=[];
            return
        end


        Moves = zeros(allotedMoves,4);
        thiscounter = 0;
        currentPos = startPos;

        while (thiscounter<allotedMoves)
            thiscounter=thiscounter+1;
            Moves(thiscounter,1:2)=currentPos;
            currentDist=sum(abs(currentPos-endPos));
            %thisMove = repmat(currentPos,[4,1])+MovingMat;
            thisMove = bsxfun(@plus,MovingMat,currentPos);
            thisPoss = ~thisWired((thisMove(:,2)-1)*size(thisWired,1)+thisMove(:,1));

            %thisDist = sum(abs(thisMove - repmat(endPos,[4,1])),2);
            thisDist = sum(abs(bsxfun(@plus,thisMove,-endPos)),2);
            
            if any(ismember(wiredInds,(currentPos(:,2)-1)*size(thisWired,1)+currentPos(:,1)))
                %thiscounter=thiscounter-1;
                break
            end
            
            if any(thisDist==0)
                Moves(thiscounter,3:4)=thisMove(thisDist==0,:);
                break
            end
            if thiscounter <= numPasses-2
                thisMask = (thisDist+~thisPoss.*10000) <= currentDist+1;
            else
                thisMask = (thisDist+~thisPoss.*10000) < currentDist;
            end

            if ~any(thisMask)
                thiscounter=0;
                break

            end
            Moves(thiscounter,3:4)=thisMove(find(thisMask,1),:);

            currentPos = Moves(thiscounter,3:4);

        end
        %if counter was reset to 0 then Moves will be empty
        Moves = Moves(1:thiscounter,:);
        
        if thiscounter==0
            Moves=FindShortestPath(Pos,order+1);
        end
        
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

        alreadyParsed=ismember(thesePairs,myConnections(:,1:5),'rows');
        if any(alreadyParsed)
            thesePairs=thesePairs(~alreadyParsed,:);
        end

    end


    function removeConnects=JunkCollector(value,strength)
        
        if nnz(myConnections(:,6:7)==value) > strength
            AllConnects=find(sum(myConnections(:,6:7)==value,2)>0);

            OtherPegs=nonzeros(myConnections(AllConnects,6:7).*(myConnections(AllConnects,6:7)~=value));

            NumConnects=sum(myID(OtherPegs,:)>0,2)+sum(myID(:,OtherPegs)>0,1)';

            toLose = NumConnects<strength;

            removeConnects=[nonzeros(myID(value,OtherPegs(toLose)));nonzeros(myID(OtherPegs(toLose),value))];

            %visualize(B,myMoves(:,1:4)-1);
            %visualize(B,myMoves(ismember(myMoves(:,5),removeConnects),1:4)-1);

            myConnections(removeConnects,:)=0;

            myMoves(ismember(myMoves(:,5),removeConnects),:)=0;
            %visualize(B,myMoves(myMoves(:,1)~=0,1:4)-1);
        end

    end
end