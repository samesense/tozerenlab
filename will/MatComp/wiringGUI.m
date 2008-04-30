function wiringGUI(B,W)
% WIRINGGUI Plays the Wiring MATLAB Contest game with a Graphic User Interface
% 
% Use the left/right-mouse button to draw/remove wires.
%
% WIRINGGUI(B)     opens the GUI with components in B and an empty solution.
% WIRINGGUI(B,W)   opens the GUI and scores the solution W.
%
% The MATLAB Contest Team 
% April, 2008 

[ro,co] = size(B);

if nargin==1 
    W = zeros(0,4);
end

hf = findall(0,'Tag','WiringGUI');
if isempty(hf)
    hf = figure('MenuBar','none','Name','The Wiring MATLAB Contest',...
        'Tag','WiringGUI','NumberTitle','off','WindowButtonUpFcn',@myUp,...
        'WindowButtonDownFcn',@myDown,'WindowButtonMotionFcn',@myMove);
else
    figure(hf);clf
end

uicontrol('Style','Text','String', 'Cost for Faults','Position', [280 30 80 15]);
uicontrol('Style','Edit','String', '0','Position', [280 10 80 20],'Tag','cmpts','Call',@updateGUI);
uicontrol('Style','Text','String', 'Copper Cost','Position', [360 30 80 15]);
uicontrol('Style','Edit','String', '0','Position', [360 10 80 20],'Tag','copper','Call',@updateGUI);
uicontrol('Style','Text','String', 'Total Cost','Position', [440 30 80 15]);
uicontrol('Style','Edit','String', '0','Position', [440 10 80 20],'Tag','cost','Call',@updateGUI);
uicontrol('String','<','Call',@myUndo,'Pos',[100 10 30 30]);
uicontrol('String','>','Call',@myRedo,'Pos',[140 10 30 30]);

S.W = W;
updateGUI
pm = []; p = [];

    function myMove(varargin)
        q = mean(get(gca,'CurrentPoint'));
        p = round(q*2)/2;
        set(findall(gcf,'Type','line','Tag','dotindicator'),'XData',p(1),'YData',p(2),'Visible','off')
        if p(1)>=1 && p(1)<= co && p(2)>=1 && p(2)<= ro && any(~rem(p([1 2]),1))
            if all(~rem(p([1 2]),1)) 
                if ~B(p(2),p(1))
                    h =any((W(:,[1 3])==p(2))&(W(:,[2 4])==p(1)),2);
                    if sum(diff(W(h,[1,3]),[],2)) && sum(diff(W(h,[2,4]),[],2))
                        set(findall(gcf,'Type','line','Tag','dotindicator'),'visible','on')
                    end
                end
            else
                set(findall(gcf,'Type','line','Tag','dotindicator'),'visible','on')
            end
        end        
    end
    function myDown(varargin)
        if strcmp(get(findall(gcf,'Type','line','Tag','dotindicator'),'visible'),'on')
            pm = p;
        else
            pm = [];
        end
    end
    function myUp(varargin)
        if ~isempty(pm)
            if all(~rem(pm([1 2]),1))
                w = pm([2 1 2 1]);
            else
                w = round(pm([2 1 2 1])+[-.1 -.1 .1 .1]);
            end
            switch get(gcf,'SelectionT')
                case 'normal'
                    if all(~all([W(:,1)==w(1) W(:,2)==w(2) W(:,3)==w(3) W(:,4)==w(4)],2)) && all(~all([W(:,1)==w(3) W(:,2)==w(4) W(:,3)==w(1) W(:,4)==w(2)],2))
                        W = [W;w];
                    end
                case 'alt'
                    if (w(1)==w(3)) && (w(2)==w(4))
                        W(((W(:,1)==w(1)) & (W(:,2)==w(2))) | ((W(:,3)==w(1)) & (W(:,4)==w(2))),:)=[];
                    else
                        W(((W(:,1)==w(1)) & (W(:,2)==w(2))) & ((W(:,3)==w(1)) & (W(:,4)==w(2))),:)=[];
                        W(((W(:,1)==w(3)) & (W(:,2)==w(4))) & ((W(:,3)==w(3)) & (W(:,4)==w(4))),:)=[];
                        W(((W(:,1)==w(1)) & (W(:,2)==w(2))) & ((W(:,3)==w(3)) & (W(:,4)==w(4))),:)=[];
                        W(((W(:,1)==w(3)) & (W(:,2)==w(4))) & ((W(:,3)==w(1)) & (W(:,4)==w(2))),:)=[];
                    end
            end
            updateGUI
            if isfield(S,'redo')
                S = rmfield(S,'redo');
            end
            S.undo = S;
            S.W = W;
        end
    end
    function myUndo(varargin)
        if isfield(S,'undo')
            S.undo.redo = S; 
            S = S.undo; 
            W = S.W;
            updateGUI
        end
    end
    function myRedo(varargin)
        if isfield(S,'redo')
            S = S.redo;
            W = S.W;
            updateGUI
        end
    end
    function updateGUI(varargin)
        [total_cost,BB] = visualize(B,W,hf);
        plot(0,0,'.r','markerSize',20,'tag','dotindicator','visible','off','linestyle','-')
        set(findall(gcf,'Type','uicontrol','Tag','cmpts'),'string',sprintf('%d',sum(BB(:))))
        set(findall(gcf,'Type','uicontrol','Tag','copper'),'string',sprintf('%d',total_cost-sum(BB(:))))
        set(findall(gcf,'Type','uicontrol','Tag','cost'),'string',sprintf('%d',total_cost))
    end
end
