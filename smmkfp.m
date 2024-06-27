function estimate = smmkfp(u,estimate,model)
%SMMKFP Predictive step of Kalman filter for switching multiple model
% estimate = SMMKFP(u,estimate,model) computes predictive mean values,
% covariance matrices, and probabilities of model sequences based on
% filtering values, model parameters, and input

% Dimension of state and number of filtering estimates
[nx,nXf] = size(estimate.xfseq);

% Number of modes (models)
nModels = length(model.M);

% Set predictive values to zeros
estimate.xpseq = zeros(nx,nXf);
estimate.Pxxpseq = zeros(nx,nx,nXf);
estimate.pmupseq = zeros(nModels*nXf,1);

% Compute predictions
idxModel = 0;
for i = 1:nXf 
    idxModel = idxModel + 1;
    if idxModel > nModels
        idxModel = 1;
    end
    
    % Compute predictive mean values and covariances for model sequences
    [estimate.xpseq(:,i),estimate.Pxxpseq(:,:,i)] = ...
        kfp(u,estimate.xfseq(:,i),estimate.Pxxfseq(:,:,i),...
        model.M(idxModel).A,model.M(idxModel).B,model.M(idxModel).G*model.M(idxModel).G',model.M(idxModel).q);
    
    % Compute the predictive probabilities of model sequences
    estimate.pmupseq((i-1)*nModels+1:i*nModels) = model.P(:,idxModel)*estimate.pmufseq(i);
end


end