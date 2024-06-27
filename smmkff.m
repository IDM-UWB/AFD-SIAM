function estimate = smmkff(y,estimate,model,computeOnlyPredictiveOutput)
%SMMKFF Filtering step of Kalman filter for switching multiple model
%   estimate = SMMKFF(y,estimate,model,computeonlypredictiveoutput)
%   computes filtering mean values, covariance matrices, and probabilities
%   of model sequences based on predictive values, model parameters,
%   and measurement

% !!!!!!!!!!!!!!!!
% Square root version should be used in future - needs to be done
% !!!!!!!!!!!!!!!!!!!!!!!!


if nargin<4
    computeOnlyPredictiveOutput = false;
end

% Dimension of state and number of predictive estimates
[nx,nXp] = size(estimate.xpseq);

% Dimension of output
ny = size(model.M(1).C,1);

% Number of modes (models)
nModels = length(model.M);

% Determined the number of filtering estimates
nXf = nModels*nXp;

% Initialize state filtering estimates by zeros
estimate.xfseq = zeros(nx,nXf);
estimate.Pxxfseq = zeros(nx,nx,nXf);
estimate.pmufseq = zeros(nXf,1);

% Initialize output predictive estimates by zeros
estimate.ypseq = zeros(ny,nXf);
estimate.Pyypseq = zeros(ny,ny,nXf);

% Initialize auxiliary variables for numerically robust computation of
% filtering probabilities of model sequences
expyp = zeros(1,nXf);
srdetPyp = zeros(1,nXf);

for i = 1:nXp
    for j = 1:nModels
        % Precompute index to shorten notation and speed up code
        idx = nModels*(i-1)+j;
        [estimate.xfseq(:,idx),estimate.Pxxfseq(:,:,idx),...
            estimate.ypseq(:,idx),estimate.Pyypseq(:,:,idx)] = ...
            kff(y,estimate.xpseq(:,i),estimate.Pxxpseq(:,:,i),model.M(j).C,model.M(j).H*model.M(j).H',model.M(j).r);
               
        % Compute auxiliary variabels used in numerically robust computation of
        % filtering probabilities of model sequences
        e = y - estimate.ypseq(:,idx);
        expyp(idx) = -0.5*e'*(estimate.Pyypseq(:,:,idx)\e);
        srdetPyp(idx) = 1/sqrt(det(estimate.Pyypseq(:,:,idx))); % Note that constant 2*pi is not included as it gets cancelled during weight computation
    end
end

if computeOnlyPredictiveOutput
    return
end

% Compute the filtering probabilities of model sequences in a numerically robust way
estimate.pmufseq = normalizeweightexp(expyp,srdetPyp.*estimate.pmupseq')';

% The weights that are less or equal to EPS(max_weight) are set to zero
% estimate.pmufseq(estimate.pmufseq<=eps(max(estimate.pmufseq))) = 0;


end