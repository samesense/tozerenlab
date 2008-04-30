function score = grade(B,W)
if size(W,2)~=4 || ~isnumeric(W) || ~isreal(W)
    error('W must be a numeric and real matrix with 4 columns')
end
[ro,co] = size(B);
if any(any(W(:,[1 3])>ro)) || any(any(W(:,[2 4])>co)) || any(any(W<1))
    error('At least one segment in W goes out of limits')
end
score = concom(B,W);
