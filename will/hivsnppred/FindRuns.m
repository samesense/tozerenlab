function [array]=FindRuns(input)
%   FindRuns
%       A helper function for MakeSNPCalls.  This function finds the number
%       of consecutive runs of ONES in a vector.  This function has been
%       vectorized so input any shape matrix.
%
%       [FORWARD_LOOK REVERSE_LOOK]=FindRuns(INPUT_ARRAY)
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
%
%
%
%


% x=input;
% [m,n] = size(input);
% reverse_looking = [zeros(1,m);input.'];
% reverse_looking = reverse_looking(:);
% p = find(~reverse_looking);
% reverse_looking(p) = [0;1-diff(p)];
% reverse_looking = reshape(cumsum(reverse_looking),[],m).';
% reverse_looking(:,1) = [];








array=zeros(size(input));
spots=find(input);
if ~isempty(spots)
    array(1:spots(1))=(spots(1)-1):-1:0;
    for k=2:length(spots)-1
        array(spots(k-1):spots(k))=(spots(k)-spots(k-1)):-1:0;
    end
    array(spots(end):end)=(length(array)-spots(end)):-1:0;

else
    array(1:end)=(length(array)-1:-1:0);
end
end
