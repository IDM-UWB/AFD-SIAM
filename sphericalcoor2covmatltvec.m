function [Pltvec,P] = sphericalcoor2covmatltvec(theta)
  %UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

  sizTheta = size(theta);

  if numel(sizTheta)~=2 
    error('cannot do')
  end

  n = (-1+sqrt(1+8*sizTheta(1)))/2;
  if fix(n)~=n || n<1 || imag(n)~=0
    error('dimension is wrong')
  end
  Pltvec = nan(sizTheta);
  for j = 1:sizTheta(2)
    L = zeros(n);

    idx = 0;
    for i = 1:n
      L(1:i,i) = spherical2cartesian(theta(idx+1:idx+i,:));
      idx = idx+i;
    end

    P = L'*L;

    Pltvec(:,j) = P(idxtril2vec([n n]));
  end

end

