function [V,dVfrobnorm_B,dVfrobnorm_T,dVfrobnorm_BT] = vvittd(problem,mV_h,V)
  %UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

  % Generate all possible inputs
  nU = cellfun(@length,problem.uGridLvl);
  % nu = length(uGridLvl);

  % The total number of all discrete inputs
  nInputs = prod(nU);

  % Generate all inputs from quantization levels given by inputs
  inputs = deaggregationidx(problem.uGridLvl,linidx2subidx(nU,1:nInputs));

  % Initialize the Bellman function
  rnkV = ones(1,length(problem.sizV)-1);
  if isempty(V)
    V = dttdzeros(problem.sizV,rnkV);
  end
  ttdCrossTol = 0.1;
  ttdCrossMaxSweep = 2;
  iter = 0;
  maxIter = 10;
  while true
    iter = iter + 1;

    getV = @(subidx) getvp(subidx',V,inputs,problem.eta,problem.mL_h,mV_h);

    z = greedy2_cross(problem.sizV, getV, ttdCrossTol,'nswp',ttdCrossMaxSweep);
    % Vnew = core2cell(z)';
    Vnew = kvadratyAFD(core2cell(z)');
    dVfrobnorm_B = chybaten(problem.sizV,V,getV,5);
    dVfrobnorm_T = chybaten(problem.sizV,Vnew,getV,5);
    dVfrobnorm_BT = dttddistfrob(V,Vnew)/sqrt(problem.gridElements);
    fprintf('VI: iter <strong>#%d</strong>, eB = %e eT = %e eBT = %e \n',iter,dVfrobnorm_B,dVfrobnorm_T,dVfrobnorm_BT);
    % IF dVfrobnorm_T>dVfrobnorm_B % Novy TTD je horsi nez puvodni TTD
    %z = dmrg_cross(length(sizV),sizV, getV, ttdCrossTol,'nswp',ttdCrossMaxSweep);
    V = Vnew;
    if iter>=maxIter
      break
    end
  end


end

