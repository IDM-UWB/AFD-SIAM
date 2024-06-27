function y = spherical2cartesian(x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sizX = size(x);

if numel(sizX)~=2
    exception = MException('MyToolbox:spherical2cartesian:wrongInputArgs',...
    'The function ''%s'' requires that the first input argument is a matrix which columns are spherical coordinates.',mfilename');
    throw(exception);
end

if any(x(1,:)<0)
    exception = MException('MyToolbox:spherical2cartesian:wrongInputArgs',...
    'The function ''%s'' requires that the first row of the input arguments is positive.',mfilename');
    throw(exception);
end

angles = x(2:sizX(1),:);
if any(angles(:)<0)|| any(angles(:)>pi)
    exception = MException('MyToolbox:spherical2cartesian:wrongInputArgs',...
    'The function ''%s'' requires that the second and other wors are in the range [0,pi].',mfilename');
    throw(exception);
end

% Preallocate array
y = nan(sizX);

tmp = x(1,:);
for i = 1:sizX(1)-1
    y(i,:) = tmp.*cos(x(i+1,:));
    tmp = tmp.*sin(x(i+1,:));
end
y(sizX(1),:) = tmp;


end