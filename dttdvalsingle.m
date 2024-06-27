function Fval = dttdvalsingle(subidx,Fttd)
%DTTDVALSINGLE evaluates TTD for a single vector of indices
%   Fval = dttdvalsingle(subidx,Fttd) returns value 'Fval' of tensor
%   represented by its TTD 'Fttd' for the column index 'idx'. The function
%   works also for tensor train decompositions that are not scalars (i.e.
%   the outer ranks r_{0} and r_{d} can be different from ones).
%
%   Example
%   F{1}(1,:,:) = [1 2;3 4; 5 6];
%   F{2} = [1 2 4 9;1 2 3 9];
%   Fval = dttdval([2;4],F)
%   Fval =
%   63
%
%   See also DTTD2FULL, DTTDFIT.
%
%   Copyright 2023 IDM (Ivo Punèocháø)


% Get sizes of input arguments
sizSubidx = size(subidx);
sizF = size(Fttd);

% Check the input arguments
if numel(sizSubidx)~=2 || sizF(2)~=sizSubidx(1) || sizSubidx(2)~=1
    exception = MException('MyToolbox:dttdval:wrongInArg',...
        'The first input argument must be a column vector of indices');
    throw(exception);
end

if ~iscell(Fttd) || numel(sizF)~=2 || sizF(1)~=1
    exception = MException('MyToolbox:dttdval:wrongInArg',...
        'The second input argument must be a row cell array');
    throw(exception);
end

% Evaluate tensor by left to right sweep
Fval = 1;
for i = 1:sizF(2)
    Fval = Fval*reshape(Fttd{i}(:,subidx(i),:),size(Fttd{i},1),size(Fttd{i},3));
end


end