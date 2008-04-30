function [score,BB] = visualize(B,W,fh)
% VISUALIZE Draws circuit components and wiring, and scores the solution
%      Inputs:
%         B - Board, a matrix with component ID's, 0's means no component.
%         W - Wiring, a matrix (nx4) matrix with wires, invalid wires are
%             ignored but charged into the score. 
%
% The MATLAB Contest Team 
% April, 2008 

if nargin==1 % No wiring
    W = zeros(0,4);
    noWiring = true;
else
    noWiring = false;
end

if nargin<3
    figure;
else
    figure(fh)
end
cla
axis equal
set(gca,'Ydir','reverse')
hold on
box on
axis([0 size(B,2) 0 size(B,1)]+.5)

[score WW Wlab C BB] = concom(B,W);

colors = hsv(max(nonzeros(C)));
m = magic(ceil(sqrt(size(colors,1))));
colors(m(m<=size(colors,1)),:)=colors;

if noWiring
   colors = colors.*0;
end
for i = 1:numel(Wlab)
    plot(WW(i,[2 4]),WW(i,[1 3]),'color',colors(Wlab(i),:))
end

[dummy,dummy,c]=unique(nonzeros(B));
colidx = B;
colidx(B>0) = c;
cmpcol = 1-(1-jet(max(c)))/4;
[i,j] = find(B);
for k = 1:numel(i)
   rectangle('Position', [j(k)-.35,i(k)-.35,.7,.7],...
       'Curvature',[.7 .7],'FaceColor',cmpcol(colidx(i(k),j(k)),:),...
       'LineWidth',1,'EdgeColor',colors(C(i(k),j(k)),:));
   if noWiring || BB(i(k),j(k))==0
       text(j(k),i(k),num2str(B(i(k),j(k))),...
           'horizontal','center');
   else
       text(j(k),i(k),num2str(B(i(k),j(k))),...
           'horizontal','center','color','r','FontWeight','bold');
   end
end










