function theta = covmatltvec2sphericalcoor(Pltvec)
  %COVMAT2SPHERICALCOOR converts SPD matrix into spherical coordinates
%   Detailed explanation goes here

  sizPltvec = size(Pltvec);

  if numel(sizPltvec)~=2 %|| sizPltvec(2)~=1
    error('incorrect input')
  end

  n = (-1+sqrt(1+8*sizPltvec(1)))/2;
  if fix(n)~=n || n<1 || imag(n)~=0
    error('dimension is wrong')
  end

  theta = nan(n,sizPltvec(2));
  for j = 1:sizPltvec(2)
    Pltvec_ = Pltvec(:,j);
    P = Pltvec_(idxvec2tril([n n]));

    sizP = size(P);
    if numel(sizP)~=2 || sizP(1)~=sizP(2) || norm(P-P')>eps
      error('The first arg should be a symmetric matrix ')
    end

    [L,flag] = chol(P);
    if flag~=0
      error('The covariance matrix must be positive definite')
    end
    idx = 0;
    for i = 1:sizP(1)
      theta(idx+1:idx+i,j) = cartesian2spherical(L(1:i,i));
      idx = idx + i;
    end
  end
end
