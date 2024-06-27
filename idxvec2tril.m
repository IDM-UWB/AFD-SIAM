function [ linidx ] = idxvec2tril(siz)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~isvector(siz) || numel(siz)>3 || numel(siz)<2
    disp('error')
end

if siz(1)~=siz(2)
    disp('error')
else
    n = siz(1);
end

if numel(siz) < 3
    m = 1;
else
    m = siz(3);
end

lengthsOfCols = n:-1:1;
cumLengthsOfCols = cumsum(lengthsOfCols);
idx = cumLengthsOfCols-lengthsOfCols+1;


linidx = zeros(siz);
for j = 1:n
    linidx(j:n,j) = idx(j):idx(j)+lengthsOfCols(j)-1;
    linidx(j,j:n) = idx(j):idx(j)+lengthsOfCols(j)-1;
end

for i = 2:m
    linidx(:,:,i) = linidx(:,:,i-1) + cumLengthsOfCols(end);
end


end

