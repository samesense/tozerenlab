function [f_look r_look]=VFindRuns(input)
%   VFindRuns
%       A helper function for MakeSNPCalls.  This function finds the number
%       of consecutive runs of ONES in a vector.  This function has been
%       vectorized so input any shape matrix.
%
%       [FORWARD_LOOK REVERSE_LOOK]=VFindRuns(INPUT_ARRAY)
%
%   INPUT = [1 1 1 0 0 0 1 0 1 0 1 1;
%            0 1 1 1 0 1 1 0 0 1 1 1];
%
%   FL =    [3 2 1 0 0 0 1 0 1 0 2 1;
%            0 3 2 1 0 2 1 0 0 3 2 1];
%
%   RL =    [1 2 3 0 0 0 1 0 1 0 1 2;
%           [0 3 2 1 0 1 2 0 0 1 2 3];
%



m = size(input,1);
 r_look = [reshape([zeros(m,1),input].',[],1);0];
 f_look = r_look;
 p = find(~r_look);
 d = 1-diff(p);
 r_look(p) = [0;d];
 r_look = reshape(cumsum(r_look(1:end-1)),[],m).';
 r_look(:,1) = [];
 f_look(p) = [d;0];
 f_look = reshape(cumsum(-f_look(1:end-1)),[],m).';
 f_look(:,end) = [];














% 
% [m,n] = size(input);
% r_look = [zeros(1,m);input.'];
% r_look = r_look(:);
% p = find(~r_look);
% r_look(p) = [0;1-diff(p)];
% r_look = reshape(cumsum(r_look),[],m).';
% r_look(:,1) = [];
% 
% f_look = [zeros(1,m);fliplr(input).'];
% f_look = f_look(:);
% p = find(~f_look);
% f_look(p) = [0;1-diff(p)];
% f_look = reshape(cumsum(f_look),[],m).';
% f_look(:,1) = [];
% f_look=fliplr(f_look);




% array=zeros(size(input));
% spots=find(input);
% if ~isempty(spots)
%     array(1:spots(1))=(spots(1)-1):-1:0;
%     for k=2:length(spots)-1
%         array(spots(k-1):spots(k))=(spots(k)-spots(k-1)):-1:0;
%     end
%     array(spots(end):end)=(length(array)-spots(end)):-1:0;
% 
% else
%     array(1:end)=(length(array)-1:-1:0);
% end