function Q = Qfcn_for_Vttd(xi,u,Vprevious_TTD,A,B,G,C,H,P,eta,gxi2theta,gtheta2subidx)
%Qfcn Q function for the active fault detection for two second order models
%with scalar measurement

% Get dimensions of output and state, and number of models
[ny,nx,nModel] = size(C);

  if nModel~=4 || ny~=1 || nx~=2
    error('not implemented')
  end

  % Select weights for the Unscented transformation
  kappa = 1;
  omega = [kappa/(1+kappa), 1/(2*(1+kappa)), 1/(2*(1+kappa))];


  % Extract statistics from xi
  for i = 1:nModel
    xf(:,i) = xi((i-1)*2+(1:2));
  end
  Pxxf(:,:,1) = [xi(9) xi(10)
  xi(10) xi(11)];
  Pxxf(:,:,2) = [xi(12) xi(13)
  xi(13) xi(14)];
  Pxxf(:,:,3) = [xi(15) xi(16)
  xi(16) xi(17)];
  Pxxf(:,:,4) = [xi(18) xi(19)
  xi(19) xi(20)];
  pmuf = [xi(21);xi(22);xi(23);1-sum(xi(21:23))];

  % Perform prediction and compute Sigma points
  xp = nan(nx,nModel);
  Pxxp = nan(nx,nx,nModel);
  yp = nan(ny,nModel,nModel);
  Pyyp = nan(ny,ny,nModel,nModel);
  Pmupf = nan(nModel,nModel);
  ySigmaPoint = nan(ny,2*ny+1,nModel,nModel);
  for i = 1:nModel
    xp(:,i) = A(:,:,i)*xf(:,i) + B(:,i)*u;
    Pxxp(:,:,i) = A(:,:,i)*Pxxf(:,:,i)*A(:,:,i)' + G(:,:,i)*G(:,:,i)';
    for j = 1:nModel
      yp(:,j,i) = C(:,:,j)*xp(:,i);
      Pyyp(:,:,j,i) = C(:,:,j)*Pxxp(:,:,i)*C(:,:,j)' + H(:,:,j)*H(:,:,j)';
      Pmupf(j,i) = P(j,i)*pmuf(i);
      tmp = sqrt(1+kappa)*chol(Pyyp(:,:,j,i),'lower');
      ySigmaPoint(:,:,j,i) = yp(:,j,i) + [zeros(ny,1) tmp -tmp];
    end
  end

  % Compute the Q function
% Q(xi,u) = L(xi) + \eta* E{ V(f(xi,u,w)) | xi,u}

  % Initialize the value of the Q function by the immediate cost L(xi)
  Q = min(pmuf);

  % Perform filtering step and compute mean of V using UT
  xfNext = nan(nx,nModel,nModel);
  PxxfNext = nan(nx,nx,nModel,nModel);
  lambda_exp = nan(nModel,nModel);
  lambda_c = nan(nModel,nModel);
  for i = 1:nModel
    for j = 1:nModel
      for ell = 1:2*ny+1
        % At this point we have particular measurement represented by
            % the sigma point ySigmaPoint(:,ell,j,i) and we perform estimation
        for iest = 1:nModel
          for jest = 1:nModel
            K = Pxxp(:,:,iest)*C(:,:,jest)'/Pyyp(:,:,jest,iest);
            e = ySigmaPoint(:,ell,j,i) - yp(:,jest,iest);
            xfNext(:,jest,iest) = xp(:,iest) + K*e;
            PxxfNext(:,:,jest,iest) = (eye(nx) - K*C(:,:,jest))*Pxxp(:,:,iest);
            lambda_exp(jest,iest) = -0.5*e/Pyyp(:,:,jest,iest)*e;
            lambda_c(jest,iest) = 1/sqrt(det(Pyyp(:,:,jest,iest)));
          end
        end
        % Pmufs = Pmupf.*lambda;
        % Pmufs = Pmufs./sum(Pmufs(:));
        
        e = lambda_exp(:)+log(Pmupf(:));
        c = lambda_c(:);
        maxE = max(e);
        % Compute non-normalized weights
        tmp = exp(e-maxE);
        % Check for possile infs, if they occur at positions where c is
        % zero it is ok an set them to zero
        flagInf = isinf(tmp);
        flagZero = c == 0;
        if all(flagZero | ~flagInf)
            tmp(flagInf) = 0;
        end

        w = c.*tmp;
        sumW = sum(w);
        % Normalize weights
        w = w./sumW;

        % Meging
        Pmufs=reshape(w,nModel,nModel);
        PmufNext = sum(Pmufs,2);

        xfNextMerge = zeros(nx,nModel);
        PxxfNextMerge = zeros(nx,nx,nModel);
        for jest = 1:nModel
          for iest = 1:nModel
            xfNextMerge(:,jest) = xfNextMerge(:,jest) + (Pmufs(jest,iest)/PmufNext(jest))*xfNext(:,jest,iest);
          end
          for iest = 1:nModel
            PxxfNextMerge(:,:,jest) = PxxfNextMerge(:,:,jest) + (Pmufs(jest,iest)/PmufNext(jest))*...
            (PxxfNext(:,:,jest,iest) + (xfNext(:,jest,iest)-xfNextMerge(:,jest))*(xfNext(:,jest,iest)-xfNextMerge(:,jest))');
          end

          [L,flag] = chol(PxxfNextMerge(:,:,jest));
          if flag~=0            
            warning('The covariance matrix must be positive definite')
            keyboard
          end
        end

        % Put statistics into vector xiNext that is the argument of the
            % Bellman function Vold
        xiNext = zeros(23,1);
        xiNext(1:2) = xfNextMerge(:,1);
        xiNext(3:4) = xfNextMerge(:,2);
        xiNext(5:6) = xfNextMerge(:,3);
        xiNext(7:8) = xfNextMerge(:,4);
        xiNext(9:10) = PxxfNextMerge(:,1,1);
        xiNext(11) = PxxfNextMerge(2,2,1);
        xiNext(12:13) = PxxfNextMerge(:,1,2);
        xiNext(14) = PxxfNextMerge(2,2,2);
        xiNext(15:16) = PxxfNextMerge(:,1,3);
        xiNext(17) = PxxfNextMerge(2,2,3);
        xiNext(18:19) = PxxfNextMerge(:,1,4);
        xiNext(20) = PxxfNextMerge(2,2,4);
        xiNext(21) = PmufNext(1);
        xiNext(22) = PmufNext(2);
        xiNext(23) = PmufNext(3);

        idxNext = gtheta2subidx(gxi2theta(xiNext));
        % idxNext
        % Continue in computation of the mean value
            % Q = Q + eta*omega(ell)*Pmupf(j,i)*Vprevious(xiNext);
        Q = Q + eta*omega(ell)*Pmupf(j,i)*dttdndlininterp(idxNext,Vprevious_TTD);
      end
    end
  end


end

