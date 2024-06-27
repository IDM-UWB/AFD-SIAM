function xiNext = phismmkf(xi,u,yNext,model,vec2tril,tril2vec,reducedProbabilityVector)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

nx = size(model.M(1).A,1);
nModels = length(model.M);

% Extract means, covariance matrices and probabilities from the hyperstate
[estimate.xfseq,estimate.Pxxfseq,estimate.pmufseq] = hyperstate2sufstat(xi,nx,nModels,vec2tril,reducedProbabilityVector);

% Compute predictions
estimateNext = smmkfp(u,estimate,model);

% Compute filtering
estimateNext = smmkff(yNext,estimateNext,model,0);

% Perform merging
[pmufseqmerged,xfseqmerged,Pxxfseqmerged] = mmmerge(estimateNext.pmufseq,estimateNext.xfseq,estimateNext.Pxxfseq,nModels);

% Put filtering means, covariance matrices and probabilities back to
% hyperstate
xiNext = sufstat2hyperstate(xfseqmerged,Pxxfseqmerged,pmufseqmerged,tril2vec,reducedProbabilityVector);


end

