function [prob] = hypercube2probvec(theta)
  %HYPERCUBE2PROBVEC converts hypercube into probability vector
%   Detailed explanation goes here

  siztheta = size(theta);

  if numel(siztheta)~=2 %|| siztheta(2)~=1
    error('incorrect input')
  end

  if any(theta<0|theta>1,'all')
    error('not a hypercube')
  end

  prob = nan(siztheta);
  for j = 1:siztheta(2)
    for i = 1:siztheta(1)
      prob(i,j) = theta(i,j)*prod(1-theta(i+1:end,j));
    end
  end
end
