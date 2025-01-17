function [pmuseqnew,idxMaxProbability] = ps2pfs(pmuseq,nSequencesNew)
%PS2PFS Computes probabilities of terminal subsequences
%   [pmuseqnew,idxMaxProbability] = ps2pfs(pmuseq,nSequencesNew) computes probabilities of
%   terminal
%   subsequences 'pmuseqterminal' given the probabilities of individual
%   sequences 'pmuseq' and the number of required terminal sequences 'n'.
%   Note that the algorithm do not check if nSequencesNew is an integer
%   power of the number of models.


% Number of model sequences
nSequences = length(pmuseq);

% Check that the number of required sequences is the same or lower than
% current number of sequnces
if nSequencesNew < nSequences
    % Reshape the column vector of probabilities to get matrix P, where
    % each column of P contains probabilities of model sequences with the
    % same final subsequence
    pmuseqReshaped = reshape(pmuseq,nSequencesNew,[]);
    
    % The probabilities of model sequences that have the same final sequence
    pmuseqnew = sum(pmuseqReshaped,2);
    
    % Index of model sequence that has maximum probability within group of sequences with common final model sequence
    [~,idxMaxProbability] = max(pmuseqReshaped,[],2);
    idxMaxProbability = (idxMaxProbability-1)*nSequencesNew + (1:nSequencesNew)';
else
    pmuseqnew = pmuseq;
    [~,idxMaxProbability] = max(pmuseq,[],2);
end


end