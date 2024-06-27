function tf = isaninteger(x)
%ISANINTEGER returns true for integer value
%   tf = isaninteger(x) returns logival array with true for integer values

tf = isfinite(x) & floor(x)==x;


end