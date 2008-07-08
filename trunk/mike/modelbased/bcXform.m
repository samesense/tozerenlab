                    				
function [XFORM,LAMBDA]=bcXform(X)      				
for i=1:size(X,1)				    				
	[XFORM(i,:),LAMBDA(i)]=boxcox(double(X(i,:))');		
end