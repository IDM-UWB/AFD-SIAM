function  linidx = idxtril2vec(siz)
%UNTITLED returns linear index of entries in the lower triangular matrix
%   linidx = idxtril2vec(siz) returns the linear index 'linidx' of elements
%   in a lower triangular marices of size 'siz'. It also works for three
%   dimensional matrices, where matrices are assumed to be stacked along
%   the third dimension.
%
%   Example
%   siz = [3 3 2];
%   idxtril2vec(siz)
%   ans =
%
%     1
%     2
%     3
%     5
%     6
%     9
%    10
%    11
%    12
%    14
%    15
%    18
%
%   See also IDXVEC2TRIL
%
%   Copyright 2015 IDM (Ivo Punèocháø)

% The code was simplified in 2023

% Get the size of the input argument
sizSiz = size(siz);

if numel(sizSiz)~=2 || (sizSiz(1)>1 && sizSiz(2)>1) || numel(siz)>3 || numel(siz)<2 || any(fix(siz)~=siz)
    exception = MException('MyToolbox:msteppred:PredTooShort',...
    'The function ''%s'' requires that the input argument is a vector of integers that has dimension 2 or 3.',mfilename');
throw(exception);
end

if siz(1)~=siz(2)
    exception = MException('MyToolbox:msteppred:PredTooShort',...
    'The function ''%s'' requires that the input argument is a vector that has the first two elements the same.',mfilename');
    throw(exception);
end

% If only 2D matrix is provided set the third dimesion to 1
if numel(siz) < 3
    siz(3) = 1;
end

% Create the linear index for the matrix in the first 'plane'
nEntriesInOneMatrix = siz(1)^2;
nEntriesInOneLowerMatrix = (siz(1)^2+siz(1))/2;

linidx = nan(siz(3)*nEntriesInOneLowerMatrix,1);
idx = 1;
for i = 1:siz(1)
     linidx(idx:idx+siz(1)-i) = ((i-1)*siz(1)+i:i*siz(1))';
     idx = idx+siz(1)-i+1;
end

% Create the linear index for the other plane if there are any
for i = 2:siz(3)
    linidx((i-1)*nEntriesInOneLowerMatrix+1:i*nEntriesInOneLowerMatrix) = linidx(1:nEntriesInOneLowerMatrix) + (i-1)*nEntriesInOneMatrix;
end

% Original code that was replaced in 2023 by much simpler code given above 
% 
% lengthsOfCols = n:-1:1;
% cumLengthsOfCols = cumsum(lengthsOfCols);
% lengthOfRC = cumLengthsOfCols(end);
% 
% idxCol = zeros(1,lengthOfRC);
% idxCol(cumLengthsOfCols-lengthsOfCols + 1) = 1;
% c = cumsum(idxCol);
% 
% idxRow = ones(1,lengthOfRC);
% idx = cumLengthsOfCols(2:end)-lengthsOfCols(2:end) + 1;
% idxRow(idx) = idxRow(idx) - lengthsOfCols(2:end);
% r = cumsum(idxRow);
% 
% linidx = zeros(lengthOfRC*m,1);
% 
% linidx(1:lengthOfRC) = sub2ind([n n],r,c);
% 
% for i = 2:m
%     linidx((i-1)*lengthOfRC+1:i*lengthOfRC) = linidx(1:lengthOfRC)+(i-1)*n^2;
% end


end

