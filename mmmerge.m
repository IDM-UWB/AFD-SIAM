function [pmuseqmerged,xseqmerged,Pxxseqmerged] = mmmerge(pmuseq,xseq,Pxxseq,nSequencesMerged)
%MMMERGE Computes statistics of terminal subsequences
%   Detailed explanation goes here

% The dimension of the state  and the number of model sequences
[nx,nSequences] = size(xseq);

if nSequencesMerged < nSequences
    % Compute probabilities of model sequences after merging
    pmuseqmerged = ps2pfs(pmuseq,nSequencesMerged);
    
    % Compute the mean values of the merged sequences
    xseqmerged = reshape(sum(reshape(bsxfun(@times,xseq,pmuseq'),nSequencesMerged*nx,[]),2),nx,nSequencesMerged);
    xseqmerged = bsxfun(@rdivide,xseqmerged,pmuseqmerged');
    
    % Replace nan that results from zero division by zeros
    xseqmerged(isnan(xseqmerged)) = 0;
    
    
    % Compute the covariance matrices of the merged sequences
    % Preallocate arays
    Pw = zeros(nx,nx,nSequences);
    Pxxseqmerged = zeros(nx,nx,nSequencesMerged);
    
    dxseq = xseq - repmat(xseqmerged,1,nSequences/nSequencesMerged);
    for i = 1:nSequences
        Pw(:,:,i) = (Pxxseq(:,:,i) + dxseq(:,i)*dxseq(:,i)')*pmuseq(i);
    end
    
    for i = 1:nSequencesMerged
        if pmuseqmerged(i) > 0
            Pxxseqmerged(:,:,i) = sum(Pw(:,:,i:nSequencesMerged:end),3)/pmuseqmerged(i);
        else
            % Prevent getting nan as a result of division by zero, set up
            % to the identity matrix. The matrix can be chosen arbitrary as
            % the probability of the sequqnce is zero anyway
            Pxxseqmerged(:,:,i) = eye(nx,nx);
        end
    end
else
    pmuseqmerged = pmuseq;
    xseqmerged = xseq;
    Pxxseqmerged = Pxxseq;
end


end