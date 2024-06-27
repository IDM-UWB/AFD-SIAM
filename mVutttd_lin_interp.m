function mVxu = mVutttd_lin_interp(subidx,u,V,varphi_h,pwxuGSParams_h,kappa)
% MVUT_TTD_LIN_INTERP  Mean value of the Bellman function when random variable has
% Gaussian sum pdf and value function is approximated using the tensor train decomposition.
%
% mVxu = MVUT_TTD_LIN_INTERP(x,u,V,inner_edges,f_h,pwxuGSParams_h,kappa) returns the
% mean value of V(f(x,u,w)) conditioned by x and u for the given state 'x'
% and input 'u'. The Bellman function V is defined over a discrete grid
% with inner edges given in 'innerEdges' and approximated using tensor
% train decomposition with tensor cores in the cell array 'V'. Function
% handles 'f_h' and 'pwxuGSParams_h' contain functions that compute new
% state and parameters of the Gaussian sum of noise w, respectivelly.
%          ╭──────────────────────────────────────────────────────────╮
%          │         uses linear interpolation/extrapolation          │
%          ╰──────────────────────────────────────────────────────────╯

% Set default value of Kappa for the unscented transformation
  if nargin < 6
    kappa = 3;
  end

  % Get parameters of Guassian sum pdf p_{w|x,u}
  [weight,mw,Pww] = pwxuGSParams_h(subidx,u);

  % Find indices of nonzero weights
  idxNonZeroWeight = find(weight>0);

  % The number of terms in GS
  nweight = length(weight);

  % Preallocate array
  Vxuw = zeros(nweight,1);

  % Compute mean value for each non-zero weight term of the Guassian sum
  for i = idxNonZeroWeight
    % Generate sigma points for the i-th term of gaussian sum
    [wSigmaPts,weighSigmaPts] = generatesigmapoints(mw(:,i),Pww(:,:,i),kappa);

    % Compute approximate mean of the i-th term of the
    for j = 1:size(wSigmaPts,2)
      xi_idx = varphi_h(subidx,u,wSigmaPts(:,j)); % necelociselny index
      V_val = 1;     %%% zde x=f(s,u), necelociselny vektor. Funguje pro rovnomerne vzorkovani.
      for ii = 1:length(xi_idx)
        % aux=(xi_idx(ii)-sNodes{ii}(1))/(sNodes{ii}(2)-sNodes{ii}(1));
        aux=xi_idx(ii)-1;
        % idx=1+floor(aux);
        %idx=floor(xi_idx(ii));
        %idx=max([idx,1]);
        %ldx=length(sNodes{ii});
        %idx=min([idx,ldx-1]);
        idx = min(max(floor(xi_idx(ii)),1),size(V{ii},2)-1);
        vaha = idx-aux;
        matice = V{ii}(:,idx,:)*vaha+V{ii}(:,idx+1,:)*(1-vaha);
        V_val = V_val*reshape(matice,size(V{ii},1),[]);   
      end
      % Vxuw(i) =  Vxuw(i) + weighSigmaPts(j)*dttdvalsingle(varphi_h(subidx,u,wSigmaPts(:,j)),V);
      Vxuw(i) =  Vxuw(i) + weighSigmaPts(j)*V_val;
    end
  end

  % Mean value E{V(f(x,u,w)|x,u}
  mVxu = weight*Vxuw;
