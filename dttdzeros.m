function Fttd = dttdzeros(sizF,rnkF)
%DTTDZEROS return tensor with zero cores
%   Fttd = dttdzeros(sizF,rnkF) returns a row cell array 'Fttd' that contains
%   TTD decomposition of a tensor for given size of tensor 'sizF' and ranks
%   'rnkF'. The elements of the tensor cores are set to zeros.
%
%   Example
%   sizF = [5 6 4];
%   rnkF = [2 3];
%   Fttd = dttdzeros(sizF,rnkF)
%   Fttd =
%   1×3 cell array
%
%    {1×5×2 double}    {2×6×3 double}    {3×4 double}
%
%   See also DTTDONES.
%
%   Copyright 2023 IDM (Ivo Punèocháø)

% Get sizes of input arguments
sizSizF = size(sizF);
sizRnkF = size(rnkF);

% Check input arguments
if numel(sizSizF)~=2 || sizSizF(1)~=1
    error('The first argument should be a row vector of sizes of the tensor.')
end

if numel(sizRnkF)~=2 || sizRnkF(1)>1 || ~( sizRnkF(2)+1==sizSizF(2) || sizRnkF(2)==sizSizF(2)+1 )
   error('The second argument should be a row vector of ranks of the tensor. The number of ranks should the number of dimensions minus one or plus one if boundary ranks are given.')
end 

% If boundary ranks are not given augment ranks by 1 at the beginnnig and end to have a single FOR loop
if sizRnkF(2)+1==sizSizF(2)
    rnkF = [1 rnkF 1]; % Default boundary ranks are ones
end
% Preallocate cell array
Fttd = cell(1,sizSizF(2));

% Generate trains of the tensor by drawing independenat samples from
% standard normal distribution
for i = 1:sizSizF(2)
    Fttd{i} = zeros([rnkF(i) sizF(i) rnkF(i+1)]);
end


end