function subidx = aggregationidx(x,xqtie)
% AGGREGATIONIDX returns subindices of grid cell for vectors
%   subidx = AGGREGATIONIDX(x,xqtie) returns subindices 'subidx' of a nD
%   grid cell with inner edges stored in the cell array 'xqtie' for vectors
%   'x' (it can be matrix of vectors given column-wise).
%
%   Example
%   xgtie = {[-1 -0.5 0 0.5 1],[-1 -0.5 0]};
%   x = [0.4;-0.2];
%   subidx = aggregationidx(x,xgtie)
%   subidx = 
%
%     4
%     3
%
%   See also DEAGGREGATIONIDX, SUBIDX2LINIDX, LINIDX2SUBIDX.
%
%   Copyright 2018 IDM (Ivo Punèocháø)


% The dimension and the number of the state vectors
sizX = size(x);
sizXqtie = size(xqtie);

if numel(sizX)~=2 || any(sizX==0)
    exception = MException('MyToolbox:aggregationidx:InconsistentDimensions',...
        'The input X must be at most a two dimensional array');
    throw(exception);
end

% Check that the dimensions of the state and inner edges are the same
if  ~iscell(xqtie) || ~((sizXqtie(1)==1 && sizXqtie(2)==sizX(1)) || (sizXqtie(2)==1 && sizXqtie(1)==sizX(1)))
    exception = MException('MyToolbox:aggregationidx:InconsistentDimensions',...
        'The dimension of the state and inner edges must be the same');
    throw(exception);
end

% Preallocate array for subscript indexes
subidx = zeros(sizX);

% Computes subscript indexes for each dimension
% NOTE: HISTC uses comparison edges(k)<= x < edge(k+1). 
% See help to HISTC for more details
% for i = 1:nx
%     [~,subidx(i,:)] = histc(x(i,:),[-inf xqtie{i} inf]);
% end

% Since HISTC is not recommended in new versions of MATLAB, function
% DISCRETIZE is used instead
for i = 1:sizX(1)
   subidx(i,:) = discretize(x(i,:),[-inf xqtie{i} inf]);
end


end