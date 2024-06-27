function [X,w,zmX] = generatesigmapoints(m,P,varargin)
%GENERATESIGMAPOINTS Generate sigma points given mean and covariance matrix
%   [X,w,zmX] = generatesigmapoints(m,P,kappa,sqrtP) generates sigma points
%   given columnwise in matrix X and corresponding weights in a row vector
%   w for given mean value M, and covariance matrix P. Optional parametrs
%   are parameter kappa (default is 3) and square root of covariance matrix
%   P (default empty).


% The dimensions input arguments
sizM = size(m);
sizP = size(P);

if numel(sizM)~=2 || sizM(2)~=1
    exception = MException('MyToolbox:generatesigmapoints:wrongInArgs',...
    'The first input agument must be a column vector');
throw(exception);
end

if numel(sizP)~=2 || sizP(1)~=sizP(2) || sizP(1)~=sizM(1)
    exception = MException('MyToolbox:generatesigmapoints:wrongInArgs',...
    'The second input argument has to be a covariance matrix');
throw(exception);
end

% Accept only two optional inputs at most
numVarArgsMax = 2;
if length(varargin) > numVarArgsMax
    error('mytoolbox:sigmapoints:TooManyInputs', ...
        'The fuction ''%s'' requires at most %i optional inputs: the scaling parameter kappa and square root of covariance matrix.', ...
        mfilename,numVarArgsMax);
end

% Set default values for optional input arguments
optArgs = {3,[]};

% Flag all nonempty input arguments in varargin
flagNotEmpty = cellfun(@isnotempty,varargin);

% Overwrite default values by the ones specified in varargin
optArgs(flagNotEmpty) = varargin(flagNotEmpty);

% Place optional args in memorable variable names
[kappa, sqrtP] = optArgs{:};

% If empty, compute sqrtP
if isempty(sqrtP)
    sqrtP = chol(P,'lower');
end

% Compute the number of sigma points
nSigmaPoints = 2*sizM(1) + 1;

% Compute the coefficient c
c = sqrt(sizM(1) + kappa);

% Create zero mean sigma points
% Variant 1 - time 0.027
zmX = [zeros(sizM(1),1) c*sqrtP -c*sqrtP];
% Variant 2 - time 0.038
%zmX = zeros(nx,2*nx+1);
%zmX(:,2:2*nx+1) = [c*sqrtPxx -c*sqrtPxx];
% Variant 3 - time 0.035
%zmX = zeros(nx,2*nx+1);
%zmX(:,2:nx+1) = c*sqrtPxx;
%zmX(:,nx+2:2*nx+1) = -c*sqrtPxx;


% Compute sigma points
X = m + zmX; % fastest, but works only in new versions of MATLAB
%X = mx(:,ones(1,nSigmaPoints)) + zmX; % second fastest, bsxfun can be
%faster if the number of sigma points is huge, but it is now expected in
%this function
%X = bsxfun(@plus,zmX,mx); % third fastest

% !!!!!!!! Due to the finite precision, the sigma points need not be
% exactly symmetrical around mean, thus (X-m)-zmX~=0

% Compute weights of sigma points
w = [kappa 0.5*ones(1,nSigmaPoints-1)]/(sizM(1) + kappa);