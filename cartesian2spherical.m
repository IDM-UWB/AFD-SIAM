function y = cartesian2spherical(x)
%CARTESIAN2SPHERICAL converts cartesian coordinate into spherical ones 
%   y = cartesian2spherical(x) computes spherical coordinates column-wise
%   in 'y' for given cartesian coordinates 'x' given column-wise. The
%   cartesian coordinates must differ from zero (origin)  and the last
%   coordinate must be positive (it is expected to be used in the
%   conjunction with the Cholesky decomposition). Each column of 'y'
%   contains the length and then angles in radians.


% Get dimensions
sizX = size(x);

% Check dimensions
if numel(sizX)~=2
    exception = MException('MyToolbox:cartesian2spherical:wrongInputArgs',...
        'The function ''%s'' requires that the first input argument is a matrix whose columns are cartesian coordinates.',mfilename');
    throw(exception);
end

% Check that the last cartesian coordinate is positive -- it is
% implementation constraint
if any(x(sizX(1),:)<0)
    exception = MException('MyToolbox:cartesian2spherical:wrongInputArgs',...
        'The function ''%s'' is implemented only to the case where the last cartesian coordinate is positive (it is assumed to be used in conjunction with Cholesky decomposition).',mfilename');
    throw(exception);
end

% Preallocate arrays
y = nan(sizX);

% Compute the norms
y(1,:) = sqrt(sum(x.^2,1));

% The origin does not have angles uniquely defined
if any(y(1,:) == 0)
    exception = MException('MyToolbox:cartesian2spherical:wrongInputArgs',...
        'The coordinates of origin cannot be uniquely converted to the spaherical coordinates.',mfilename');
    throw(exception);
end

% Normalize the coordinate
x = x./y(1,:);

% Compute angles
for j = 1:sizX(1)-1
    % Compute angles
    y(j+1,:) = acos(x(j,:));
    
    % Update remaining coordinates
    x(j+1:sizX(1),:) = x(j+1:sizX(1),:)./sin(y(j+1,:));
end

if abs(x(j+1:sizX(1),:)-1)>sqrt(eps)
    exception = MException('MyToolbox:cartesian2spherical:outOfTolerance',...
        'The accuracy of computation is out of tolerance.',mfilename');
    throw(exception);
end



end

