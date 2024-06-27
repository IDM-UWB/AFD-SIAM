function x = deaggregationidx(xqtlvl,subidx)
%DEAGGREGATIONIDX returns vectors from quantization levels
%   x = deaggregationidx(xqtlvl,subidx) creates 2D array 'x' of column vectors using
%   quantization levels in row cell array XQTLVL and indices SUBIDX
%
%   Example
%   xgtlvl = {[-1 -0.5 0 0.5 1],[-1 -0.5 0]};
%   subidx = [1;2];
%   x = deaggregationidx(xgtlvl,subidx)
%   x = 
%
%     -1
%     -0.5
%
%   See also AGGREGATIONIDX, SUBIDX2LINIDX, LINIDX2SUBIDX.
%
%   Copyright 2018 IDM (Ivo Punèocháø)

% Get dimensions of input arguments
sizSubidx = size(subidx);
sizXqtlvl = size(xqtlvl);

if numel(sizSubidx)~=2 || any(~isaninteger(subidx(:))) || any(subidx(:)<=0)
    exception = MException('MyToolbox:deaggregationidx:inconsistentDimension',...
        'The SUBIDX must be at most two dimensional array.');
    throw(exception);
end

% Check that dimensions are consistent
if ~iscell(xqtlvl) || ~((sizXqtlvl(1)==1 && sizXqtlvl(2)==sizSubidx(1)) || (sizXqtlvl(2)==1 && sizXqtlvl(1)==sizSubidx(1)))
    exception = MException('MyToolbox:deaggregationidx:inconsistentDimension',...
        'The dimensions of state and subscript must be the same.');
    throw(exception);
end

% Preallocate array
x = zeros(sizSubidx);

% Create states using subscripts
for i = 1:sizSubidx(1)
    x(i,:) = xqtlvl{i}(subidx(i,:));
end


end
