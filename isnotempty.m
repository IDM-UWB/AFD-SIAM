function tf = isnotempty(x)
%ISNOTEMPTY True for non empty array.
%   tf = ISNOTEMPTY(x) returns 1 if x is a non empty array and 0 otherwise.
 
tf = numel(x)~=0;


end