function xqtie = generateinneredges(xqtlvl)
%GENERATEINNEREDGES returns symmetric inner edges for given grid points
%   xqtie = generateinneredges(xqtlvl) returns cell array 'xqtie' of the
%   same size as 'xqtlvl' that contains the inner edges
%
%   Example
%   xqtlvl = {[-1 -0.5 0 0.5 1],[-1 -0.5 0]};
%   xqtie = generateinneredges(xqtlvl);
%   xqtie{1} 
%   ans =
%   -0.75       -0.25       0.25        0.75    
%   xqtie{2}
%   ans =
%   -0.75       -0.25
%
%   See also DEAGGREGATIONIDX, SUBIDX2LINIDX, LINIDX2SUBIDX.
%
%   Copyright 2018 IDM (Ivo Punèocháø)

% Get sizes
sizXqtlvl = size(xqtlvl);

if ~iscell(xqtlvl) || numel(sizXqtlvl)~=2 || (sizXqtlvl(1)>1 && sizXqtlvl(2)>1)
    exception = MException('MyToolbox:generateinneredges:wrongDimInArg',...
        'The function ''%s'' requires that its arguments is a vector (column or row) cell array.',mfilename');
    throw(exception);
end


% Preallocate cell array
xqtie = cell(sizXqtlvl);

% Compute inner edges for each element of the vector state
for i = 1:max(sizXqtlvl)
    xqtie{i} = xqtlvl{i}(1:end-1) +...
        0.5*(xqtlvl{i}(2:end) - xqtlvl{i}(1:end-1));
end


end