function [V,dVfrobnorm,dVfrobnorm_ci,dVfrobnorm_li] = vvittd_11D_example(thetaGridLvl,thetaGridInnerEdge)
  %UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

  % Generate all possible inputs
  % nU = cellfun(@length,uGridLvl);
  % nu = length(uGridLvl);

  % The total number of all discrete inputs
  % nInputs = prod(nU);

  % Generate all inputs from quantization levels given by inputs
  % inputs = deaggregationidx(uGridLvl,linidx2subidx(nU,1:nInputs));
 gridElements = prod(cellfun(@length,thetaGridLvl));
  
  A(:,:,1) = [0.9884,0.0458; -0.4492,0.8327];
  A(:,:,2) = [0.9884,0.0462; -0.4536,0.8496];
  B(:,:,1) = [0.0006; 0.0231];
  B(:,:,2) = [0.0006; 0.0231];
  G(:,:,1) = 0.001*eye(2);
  G(:,:,2) = 0.001*eye(2);
  C(:,:,1) = [1 0];
  C(:,:,2) = [1 0];
  H(:,:,1) = 0.001;
  H(:,:,2) = 0.001;
  P = [0.98 0.01; 0.02 0.99];

  eta = 0.9;

  % Define step of possible input values - implementation assumed exactly
  % three values
  uSet = [-10 0 10];

  % Initialize the Bellman function
  sizV = cellfun(@length,thetaGridLvl);
  rnkV = ones(1,length(sizV)-1);
  V = dttdzeros(sizV,rnkV);

  gsubidx2theta = @(subidx) deaggregationidx(thetaGridLvl,subidx);
  gtheta2xi = @(theta) [theta(1:2)
    theta(1:2)-theta(3:4)
    sphericalcoor2covmatltvec(theta(5:7))
    sphericalcoor2covmatltvec(theta(8:10))
    theta(11)];
  gxi2theta = @(xi) [xi(1:2)
    xi(1:2)-xi(3:4)
    covmatltvec2sphericalcoor(xi(5:7))
    covmatltvec2sphericalcoor(xi(8:10))
    xi(11)];
  % constant interpolation
  gtheta2subidx_ci = @(theta) aggregationidx(theta,thetaGridInnerEdge); % picks nearest neighbor
  % gtheta2subidx = @(theta) aggregationidx(theta,thetaGridInnerEdge); % picks nearest neighbor
  % gtheta2subidx = @(theta) aggregationidx_reg_grid(theta,thetaGridInnerEdge); % works only for regular grids
  % linear interpolation
  gtheta2subidx_li = @(theta) theta2realidx(theta,thetaGridInnerEdge); % real idx

  ttdCrossTol = 0.01;
  ttdCrossMaxSweep = 200;
  iter = 0;
  maxIter = 10;
  while true
    iter = iter + 1;

    % constant interpolation
    % Define Q function for particular parameters
    Qtmp1fcn_for_Vttd_ci = @(subidx,u,Vold_TTD) Qfcn_for_Vttd(gtheta2xi(gsubidx2theta(subidx)),u,Vold_TTD,A,B,G,C,H,P,eta,gxi2theta,gtheta2subidx_ci);
    Qtmp2fcn_for_Vttd_ci = @(subidx,u) Qtmp1fcn_for_Vttd_ci(subidx,u,V);
    getV_ci = @(subidx) min([Qtmp2fcn_for_Vttd_ci(subidx',uSet(1)) Qtmp2fcn_for_Vttd_ci(subidx',uSet(2)) Qtmp2fcn_for_Vttd_ci(subidx',uSet(3))]);
    % getV = @(subidx) Qtmp2fcn_for_Vttd(subidx',uSet(1));
    % getV = @(subidx) getvp(subidx',V,inputs,eta,mL_h,mV_h);
    % linear interpolation
    % Define Q function for particular parameters
    Qtmp1fcn_for_Vttd_li = @(subidx,u,Vold_TTD) Qfcn_for_Vttd(gtheta2xi(gsubidx2theta(subidx)),u,Vold_TTD,A,B,G,C,H,P,eta,gxi2theta,gtheta2subidx_li);
    Qtmp2fcn_for_Vttd_li = @(subidx,u) Qtmp1fcn_for_Vttd_li(subidx,u,V);
    getV_li = @(subidx) min([Qtmp2fcn_for_Vttd_li(subidx',uSet(1)) Qtmp2fcn_for_Vttd_li(subidx',uSet(2)) Qtmp2fcn_for_Vttd_li(subidx',uSet(3))]);
    % getV = @(subidx) Qtmp2fcn_for_Vttd(subidx',uSet(1));
    % getV = @(subidx) getvp(subidx',V,inputs,eta,mL_h,mV_h);


    fprintf('<strong>Constant interpolation</strong>\n')
    z_ci = greedy2_cross(sizV, getV_ci, ttdCrossTol,'nswp',ttdCrossMaxSweep);
    Vnew_ci = core2cell(z_ci)';
    dVfrobnorm_ci(iter) = dttddistfrob(V,Vnew_ci);
    fprintf('Value iteration: iteration <strong>#%d</strong>, || V(i+1)-V(i) || = %e \n',iter,dVfrobnorm_ci(iter)/gridElements)

    fprintf('<strong>Linear interpolation</strong>\n')
    z_li = greedy2_cross(sizV, getV_li, ttdCrossTol,'nswp',ttdCrossMaxSweep);
    %z = dmrg_cross(length(sizV),sizV, getV, ttdCrossTol,'nswp',ttdCrossMaxSweep);
    Vnew_li = core2cell(z_li)';
    dVfrobnorm_li(iter) = dttddistfrob(V,Vnew_li);
    fprintf('Value iteration: iteration <strong>#%d</strong>, || V(i+1)-V(i) || = %e \n',iter,dVfrobnorm_li(iter)/gridElements)

    V = Vnew_li;
    dVfrobnorm = dVfrobnorm_li;

    if iter>=maxIter
      break
    end
  end


end

