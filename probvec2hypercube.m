function [theta] = probvec2hypercube(prob)
  %PROBVEC2HYPERCUBE converts probability vector into hypercube
%   Detailed explanation goes here

  sizprob = size(prob);

  if numel(sizprob)~=2 %|| sizprob(2)~=1
    error('incorrect input')
  end

  if any(prob<0,'all')
    error('negative probabilities')
  end

  if any(sum(prob,1)>1)
    error('probabilities sum larger than one')
  end

  if any(prob==0 | prob==1,'all')
    error('not a bijective function')
  end

  theta = nan(sizprob);
  for j = 1:sizprob(2)
    for i = 1:sizprob(1)
      theta(i,j) = min(prob(i,j)/(1-sum(prob(i+1:end,j))),1);
    end
  end
end
