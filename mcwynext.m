function [w,yp,Pyyp] = mcwynext(xi,u,model,vec2tril,reducedProbabilityVector)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nx = size(model.M(1).A,1);
nModels = length(model.M);

% Extract means, covariance matrices and probabilities from the hyperstate
[estimate.xfseq,estimate.Pxxfseq,estimate.pmufseq] = hyperstate2sufstat(xi,nx,nModels,vec2tril,reducedProbabilityVector);

% Compute predictions
estimateNext = smmkfp(u,estimate,model);

% Simulate realization of measurement y
estimateNext = smmkff(0,estimateNext,model,true);

w = estimateNext.pmupseq';
Pyyp = estimateNext.Pyypseq;
yp =estimateNext.ypseq;


end

