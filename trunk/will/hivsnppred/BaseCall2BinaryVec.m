function [BINARY_VEC]=BaseCall2BinaryVec(BASE_CALL_VEC)
%   BaseCall2BinaryVec
%       Converts a vector in the form of BASE_CALL to a binary expansion
%       vector for the input of clustering algorithms.
%
%   BINARY_VEC=BaseCall2BinaryVec(BASE_CALL_VEC)
%
%
%   BASE_CALL_VEC must be a 1xN vector of values between 0 and 4 (nt2int
%           format)
%
%


if size(BASE_CALL_VEC,1)~=1
    error('BaseCall2BinaryVec:NOT_ROW','BASE_CALL_VEC must be 1xN vector')
elseif max(BASE_CALL_VEC)>4||min(BASE_CALL_VEC)<0
    error('BaseCall2BinaryVec:BAD_NUMS','Values of BASE_CALL_VEC must be  between 0 and 4')
else
    map=[0 0 0 0;1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    map_fun=@(x)map(x+1,:);
    BINARY_VEC=cell2mat(arrayfun(map_fun,BASE_CALL_VEC,'uniformoutput',false));
end
end