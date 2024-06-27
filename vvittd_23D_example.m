function [V,dVfrobnorm1,dVfrobnorm2,dVfrobnorm3] = vvittd_23D_example(model,problem)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Generate all possible inputs
% nU = cellfun(@length,uGridLvl);
% nu = length(uGridLvl);

% The total number of all discrete inputs
% nInputs = prod(nU);

  % Generate all inputs from quantization levels given by inputs
% inputs = deaggregationidx(uGridLvl,linidx2subidx(nU,1:nInputs));
  nModels = length(model.M);
  for i = 1:nModels
    A(:,:,i) = model.M(i).A;
    B(:,:,i) = model.M(i).B;
    G(:,:,i) = model.M(i).G;
    C(:,:,i) = model.M(i).C;
    H(:,:,i) = model.M(i).H;
  end

  P = model.P;

  eta = problem.eta;

  % Define step of possible input values - implementation assumed exactly
% three values
  uSet = [-10 0 10];

  % Initialize the Bellman function
  sizV = problem.sizV;
  gridElements = prod(sizV);
  rnkV = ones(1,length(sizV)-1);
  V = dttdzeros(sizV,rnkV);

  gsubidx2theta = problem.gsubidx2theta;
  gtheta2xi = problem.gtheta2xi;
  gxi2theta = problem.gxi2theta; 
  % constant interpolation
% gtheta2subidx_ci = @(theta) aggregationidx(theta,thetaGridInnerEdge); % picks nearest neighbor
% gtheta2subidx = @(theta) aggregationidx(theta,thetaGridInnerEdge); % picks nearest neighbor
% gtheta2subidx = @(theta) aggregationidx_reg_grid(theta,thetaGridInnerEdge); % works only for regular grids
% linear interpolation
  gtheta2subidx = problem.gtheta2subidx; % real idx

  ttdCrossTol = 0.01;
  ttdCrossMaxSweep = 200;
  iter = 0;
  maxIter = 30;
  while true
    iter = iter + 1;

    % linear interpolation
  % Define Q function for particular parameters
    Qtmp1fcn_for_Vttd = @(subidx,u,Vold_TTD) Qfcn_for_Vttd(gtheta2xi(gsubidx2theta(subidx)),u,Vold_TTD,A,B,G,C,H,P,eta,gxi2theta,gtheta2subidx);
    Qtmp2fcn_for_Vttd = @(subidx,u) Qtmp1fcn_for_Vttd(subidx,u,V);
    getV = @(subidx) min([Qtmp2fcn_for_Vttd(subidx',uSet(1)) Qtmp2fcn_for_Vttd(subidx',uSet(2)) Qtmp2fcn_for_Vttd(subidx',uSet(3))]);
    % getV = @(subidx) Qtmp2fcn_for_Vttd(subidx',uSet(1));
  % getV = @(subidx) getvp(subidx',V,inputs,eta,mL_h,mV_h);

    RECOMPUTE_TTD = true;
    nRecomputations = 0;
    maxRecomputations = 1;
    while RECOMPUTE_TTD
      nRecomputations = nRecomputations + 1;
      z = greedy2_cross(sizV, getV, ttdCrossTol,'nswp',ttdCrossMaxSweep);
      % Vnew{nRecomputations} = core2cell(z)';
      Vnew{nRecomputations} = kvadratyAFD(core2cell(z)');
      E1(nRecomputations) = chybaten(sizV,Vnew{nRecomputations},getV,5);
      E2(nRecomputations) = chybaten(sizV,V,getV,5);
      E3(nRecomputations) = dttddistfrob(V,Vnew{nRecomputations})/sqrt(gridElements);
      fprintf('VI: iter <strong>#%d</strong>, rec <strong>#%d</strong>, E1 = %e E2 = %e E3 = %e \n',iter,nRecomputations,E1(nRecomputations),E2(nRecomputations),E3(nRecomputations));
      % IF E1>E2 % Novy TTD je horsi nez puvodni TTD
  % maximalne 5x greedy2_cross a vybrat podle E1 ten nejlepsi pokud 5x neuspeje
      % if (E1(nRecomputations) < E2(nRecomputations)) || (nRecomputations == maxRecomputations)
      if  (nRecomputations == maxRecomputations)
        RECOMPUTE_TTD = false;
      end
    end
    if nRecomputations < maxRecomputations
      idx = nRecomputations;
    else
      [~,idx] = min(E1);
    end
    strE = {'E1';'E2';'E3'};
    mE = [mean(E1);mean(E2);mean(E3)];
    vE = [var(E1);var(E2);var(E3)];
    T = table(strE,mE,vE)

    dVfrobnorm1(iter) = E1(idx); 
    dVfrobnorm2(iter) = E2(idx);
    dVfrobnorm3(iter) = E3(idx);
    clear E1 E2 E3
    V = Vnew{idx};
    if iter>=maxIter
      break
    end
  end
end
