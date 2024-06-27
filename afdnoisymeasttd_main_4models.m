close all
clear all
clc

% Add path to functions needed in this example
addpath('../dttd')

scenario = 'second-order-oscillating';
valueFunctionFileName = 'bellman23D_coarsegrid';
% simulationResultsFileName = 'simulationResults_coarsegrid';
simulationResultsFileName = 'simulationResults_passive';

switch scenario

case 'second-order-oscillating'
  % Select sampling period
  Ts = 0.05;
  m_p = 2;
  g_p = 9.81;
  ell_p = 1;
  b_p = 6;
%          ╭──────────────────────────────────────────────────────────╮
%          │                      Define model 1                      │
%          ╰──────────────────────────────────────────────────────────╯
        % continuous model 
  A = [0 1
  -g_p/ell_p -b_p/(m_p*ell_p^2)];
  B = [0;1/(m_p*ell_p^2)];
  C = [1 0];
  D = 0;
  % discretization
  [model.M(1).A,model.M(1).B,model.M(1).C] = ssdata(c2d(ss(A,B,C,D),Ts));
  model.M(1).G = 0.001*eye(2);
  model.M(1).q = zeros(2,1);
  model.M(1).H = 0.001;
  model.M(1).r = 0;

%          ╭──────────────────────────────────────────────────────────╮
%          │                      Define model 2                      │
%          ╰──────────────────────────────────────────────────────────╯
        % continuous model 
  model.M(2) = model.M(1);
  b_p = 6.8;
  A = [0 1
  -g_p/ell_p -b_p/(m_p*ell_p^2)];
  % discretization
  model.M(2).A = ssdata(c2d(ss(A,B,C,D),Ts));

%          ╭──────────────────────────────────────────────────────────╮
%          │                      Define model 3                      │
%          ╰──────────────────────────────────────────────────────────╯
        % continuous model 
  model.M(3) = model.M(1);
  b_p = 5.2;
  A = [0 1
  -g_p/ell_p -b_p/(m_p*ell_p^2)];
  % discretization
  model.M(3).A = ssdata(c2d(ss(A,B,C,D),Ts));


%          ╭──────────────────────────────────────────────────────────╮
%          │                      Define model 4                      │
%          ╰──────────────────────────────────────────────────────────╯
        % continuous model 
  model.M(4) = model.M(1);
  b_p = 7.6;
  A = [0 1
  -g_p/ell_p -b_p/(m_p*ell_p^2)];
  % discretization
  model.M(1).A = ssdata(c2d(ss(A,B,C,D),Ts));

  % Transition probability matrix
  model.P = [0.97 0.01 0.01 0.01;
  0.01 0.95 0.01 0.02;
  0.01 0.02 0.97 0.01;
  0.01 0.02 0.01 0.96];


  % The predictive mean of the state
  model.meanx0p = [0 0]';

  % The predictive covariance of the state
  model.Pxx0p = diag([7e-5 2e-4]);

  % Predictive probabilities of models
  model.pmu0p = [1-3e-3 1e-3 1e-3 1e-3]';

  % The number of models
  nModel = length(model.M);

  % The dimension of state and input
  [nx,nu] = size(model.M(1).B);

  % The dimension of measurement
  ny = size(model.M(1).C,1);

  % The dimension of state noise
  nw = size(model.M(1).G,2);

  % The dimension of measurement noise
  nv = size(model.M(1).H,2);


  % Detection cost function - not used
  problem.Ld = ones(nModel) - eye(nModel);

  % Discount factor
  problem.eta = 0.9;
  %          ╭──────────────────────────────────────────────────────────╮
%          │     Define quantization levels for intermediate grid     │
%          ╰──────────────────────────────────────────────────────────╯

  FineGrid = false;
  if FineGrid
    nLvlMean = 61;
    nLvlvar = 10;
    nLvlProb = 1e2;

    xGridLvl1 = linspace(-1,1,nLvlMean); %5
    xGridLvl2 = linspace(-2,2,nLvlMean); %21
    xGridLvl3 = linspace(-0.5,0.5,nLvlMean); %11
    xGridLvl4 = linspace(-0.5,0.5,nLvlMean); %11
    ellGridLvl11 = sqrt(linspace(0.4e-5,0.6e-5,nLvlvar)); %4
    ellGridLvl22 = sqrt(linspace(2e-5,3e-5,nLvlvar)); %4
    anglesGridLvl = linspace(0,pi,10);
    % anglesGridLvl = [pi/2-0.1 pi/2 pi/2+0.1];
    probGridLvl = linspace(1e-3,1-1e-3,nLvlProb);
  else
    xGridLvl1 = -1:0.5:1;
    xGridLvl2 = -2:0.2:2;
    xGridLvl3 = -0.5:0.1:0.5;
    xGridLvl4 = -0.5:0.1:0.5;
    ellGridLvl11 = sqrt(linspace(0.4e-5,0.6e-5,4));
    ellGridLvl22 = sqrt(linspace(2e-5,3e-5,4));
    anglesGridLvl = [pi/2-0.1 pi/2 pi/2+0.1];
    probGridLvl = linspace(1e-3,1-1e-3,100);
  end
  thetaGridLvl = {xGridLvl1,xGridLvl2,xGridLvl3,xGridLvl4,xGridLvl3,xGridLvl4,xGridLvl3,xGridLvl4,...
  ellGridLvl11,...
  ellGridLvl22,anglesGridLvl,...
  ellGridLvl11,...
  ellGridLvl22,anglesGridLvl,...
  ellGridLvl11,...
  ellGridLvl22,anglesGridLvl,...
  ellGridLvl11,...
  ellGridLvl22,anglesGridLvl,...
  probGridLvl,probGridLvl,probGridLvl};
  thetaGridInnerEdge = generateinneredges(thetaGridLvl);

  uGridLvl = {-10:10:10};

  problem.sizV = cellfun(@length,thetaGridLvl);
  %          ╭──────────────────────────────────────────────────────────╮
%          │                  Grid mapping functions                  │
%          ╰──────────────────────────────────────────────────────────╯
  problem.gsubidx2theta = @(subidx) deaggregationidx(thetaGridLvl,subidx);
  problem.gtheta2xi = @(theta) [theta(1:2,:)
  theta(1:2,:)-theta(3:4,:)
  theta(1:2,:)-theta(5:6,:)
  theta(1:2,:)-theta(7:8,:)
  sphericalcoor2covmatltvec(theta(9:11,:))
  sphericalcoor2covmatltvec(theta(12:14,:))
  sphericalcoor2covmatltvec(theta(15:17,:))
  sphericalcoor2covmatltvec(theta(18:20,:))
  hypercube2probvec(theta(21:23,:))];

  %          ╭──────────────────────────────────────────────────────────╮
%          │          Inverses of the grid mapping functions          │
%          ╰──────────────────────────────────────────────────────────╯
  problem.gxi2theta = @(xi) [xi(1:2,:)
  xi(1:2,:)-xi(3:4,:)
  xi(1:2,:)-xi(5:6,:)
  xi(1:2,:)-xi(7:8,:)
  covmatltvec2sphericalcoor(xi(9:11,:))
  covmatltvec2sphericalcoor(xi(12:14,:))
  covmatltvec2sphericalcoor(xi(15:17,:))
  covmatltvec2sphericalcoor(xi(18:20,:))
  probvec2hypercube(xi(21:23,:))];

  problem.gtheta2subidx = @(theta) theta2realidx(theta,thetaGridInnerEdge); % real idx
otherwise
  error('Unknown scenario')
end

% Auxiliary indexing matrices
tril2vec = idxtril2vec([nx,nx,nModel]);
vec2tril = idxvec2tril([nx,nx,nModel]);
reducedProbabilityVector = true;

% Compute number of states of the full state information model
if reducedProbabilityVector
  nxi = nModel*(nx+nx*(nx+1)/2+1)-1;
else
  nxi = nModel*(nx+nx*(nx+1)/2+1);
end

% Define function g2
gtheta2subidx_control = @(theta) aggregationidx(theta,thetaGridInnerEdge); % picks nearest neighbor
% gtheta2subidx = @(theta) aggregationidx(theta,thetaGridInnerEdge); % picks nearest neighbor
% gtheta2subidx = @(theta) aggregationidx_reg_grid(theta,thetaGridInnerEdge); % works only for regular grids
problem.gtheta2subidx = @(theta) theta2realidx(theta,thetaGridInnerEdge); % real idx

% Define system wraped into grid mappings
fsux_h = @(subidx,u,y) problem.gtheta2subidx(problem.gxi2theta(phismmkf(problem.gtheta2xi(problem.gsubidx2theta(subidx)),u,y,model,vec2tril,tril2vec,reducedProbabilityVector)));

% Define predicitve distribution of future output
pxsugs_h = @(subidx,u) mcwynext(problem.gtheta2xi(problem.gsubidx2theta(subidx)),u,model,vec2tril,reducedProbabilityVector);

% Define cost functions
% mL_h = @(subidx,u)  min(probGridLvl(subidx(end)),1-probGridLvl(subidx(end)));
mL_h = @(subidx,u)  min([probGridLvl(subidx(end-2:end)) 1-sum(probGridLvl(subidx(end-2:end)))]);
%mV_h = @(subidx,u,V) mVutttd(subidx,u,V,fsux_h,pxsugs_h);
mV_h = @(subidx,u,V) mVutttd_lin_interp(subidx,u,V,fsux_h,pxsugs_h);

% trainBellman = 'fresh';
%trainBellman = 'continue';
% trainBellman = 'no';
% trainBellman = 'fresh_23D_example';
% trainBellman = 'passive';
trainBellman = 'NN';
switch trainBellman
case 'fresh'
  V = vvittd(thetaGridLvl,uGridLvl,eta,mL_h,mV_h,[]);
  save('bellman','V')
case 'continue'
  load('bellman','V')
  V = vvittd(thetaGridLvl,uGridLvl,eta,mL_h,mV_h,V);
  save('bellman','V')
case 'fresh_23D_example'
  % [V,dVfrobnorm,dVfrobnorm_ci,dVfrobnorm_li] = vvittd_11D_example(thetaGridLvl,thetaGridInnerEdge);
        % [V,dVfrobnorm,dVfrobnorm_ci,dVfrobnorm_li] = vvittd_11D_example_ci(thetaGridLvl,thetaGridInnerEdge);
        % save('bellman31_ci','dVfrobnorm_li','dVfrobnorm_ci','thetaGridLvl')
  gridElements = prod(cellfun(@length,thetaGridLvl));
  [V,dVfrobnorm1,dVfrobnorm2,dVfrobnorm3] = vvittd_23D_example(model,problem);
  save(valueFunctionFileName,'V','dVfrobnorm1','dVfrobnorm2','dVfrobnorm3','gridElements','thetaGridLvl')
case 'NN'
  % load('NN2_parameters_continuous');
  problem_tmp = problem;
  load('NN2_parameters_discrete');
  problem = problem_tmp;
  x_dim = model.x_dim;
  u_dim = 1;
  neurons = policy.approx.neurons;
  W = reshape(theta(1:(x_dim+1)*neurons),(x_dim+1),neurons);
  C = reshape(theta((x_dim+1)*neurons+1:end),neurons,u_dim);
  % for continuous bounded policy
  % policy_NN = @(X) max(model.min_u,min(model.max_u,mlp_out(C,W,X,model,policy)));
  % for discrete policy
  policy_NN = @(X) round(max(model.min_u,min(model.max_u,mlp_out(C,W,X,model,policy)))/model.max_u)*model.max_u;
case 'passive'
otherwise
  load(valueFunctionFileName)
end
% figure,plot(dVfrobnorm)

%          ╭──────────────────────────────────────────────────────────╮
%          │                   input specification                    │
%          ╰──────────────────────────────────────────────────────────╯
% Generate all possible inputs
nU = cellfun(@length,uGridLvl);
nu = length(uGridLvl);

% The total number of all discrete inputs
nInputs = prod(nU);

% Generate all inputs from quantization levels given by inputs
inputs = deaggregationidx(uGridLvl,linidx2subidx(nU,1:nInputs));

controller = @(xi) getvp(gtheta2subidx_control(problem.gxi2theta(xi)),V,inputs,problem.eta,mL_h,mV_h);
%%

%          ╭──────────────────────────────────────────────────────────╮
%          │                    system simulation                     │
%          ╰──────────────────────────────────────────────────────────╯

simulateSystem = true;
nMC = 1e3; % Monte Carlo simulations
J = zeros(1,nMC);
if simulateSystem
  for iMC = 1:nMC
    % Simulate the trajectories
    F = 500;
    Fp1 = F + 1;

    % Preallocate arrays
    x = nan(nx,Fp1);
    mu = nan(1,Fp1);
    y = nan(ny,Fp1);
    d = nan(1,Fp1);
    u = nan(nu,Fp1);

    str.xpseq = zeros(nx,nModel);
    str.Pxxpseq = zeros(nx,nx,nModel);
    str.pmupseq = zeros(nModel,1);
    str.ypseq = zeros(ny,nModel);
    str.Pyypseq = zeros(ny,ny,nModel);
    str.xfseq = zeros(nx,nModel);
    str.Pxxfseq = zeros(nx,nx,nModel);
    str.pmufseq = zeros(nModel^2,1);
    estimate = repmat(str,1,Fp1);


    % Initialize estimator
    estimate(1).xpseq = model.meanx0;
    estimate(1).Pxxpseq = model.Pxx0;
    estimate(1).pmupseq = model.pmu0;

    % generate initial state x and mu
    x(:,1) = model.meanx0 + chol(model.Pxx0,'lower')*randn(nx,1);
    mu(1) = gendrnd(model.pmu0);
    % generate noises
    w = randn(nw,F);
    v = randn(nv,Fp1);
    for k = 1:Fp1
      % generate measurments
      y(:,k) = model.M(mu(k)).C*x(:,k) + model.M(mu(k)).r + model.M(mu(k)).H*v(:,k);

      % Kalman filter measurement update
      estimate(k) = smmkff(y(:,k),estimate(k),model,0);

      % merging hypotheses
      [estimate(k).pmufseq,estimate(k).xfseq,estimate(k).Pxxfseq] = mmmerge(estimate(k).pmufseq,estimate(k).xfseq,estimate(k).Pxxfseq,nModel);

      [~,d(k)]=max(estimate(k).pmufseq);


      %%%%%%%%%%%%%%%%%%%%%%%% AFD 
      switch trainBellman
      case {'no','fresh', 'continue', 'fresh_23D_example'}
        % form xi
        xi = sufstat2hyperstate(estimate(k).xfseq,estimate(k).Pxxfseq,estimate(k).pmufseq,tril2vec,reducedProbabilityVector);
        % Generate input
        [~,u(:,k)] = controller(xi);
      case 'passive'
        %%%%%%%%%%%%%%%%%%%%%%%% PFD 
        u(:,k) = 10*randi([-1 1],1,1); % for passive FD
      case 'NN'
        %%%%%%%%%%%%%%%%%%%%%%%% NN AFD
        xi = sufstat2hyperstate(estimate(k).xfseq,estimate(k).Pxxfseq,estimate(k).pmufseq,model.tril2vec,model.reducedProbabilityVector);
        u(:,k) = policy_NN(xi);
      end


      if k<Fp1
        % generate state x and mu at next time instant
        x(:,k+1) = model.M(mu(k)).A*x(:,k) + model.M(mu(k)).B*u(:,k) + model.M(mu(k)).q + model.M(mu(k)).G*w(:,k);
        mu(k+1) = gendrnd(model.P(:,mu(k)));

        % Kalman filter time update
        estimate(k+1) = smmkfp(u(:,k),estimate(k),model);
      end
    end
    % calculate criterion
    t = 0:F;
    J(iMC) = sum(problem.eta.^t.*sum(d~=mu,1));
    fprintf('%u/%u MC simulations done.\n',iMC,nMC);
  end
  % save('simulationResults_finegrid_P1st','J');
  save(simulationResultsFileName,'J');
  % AFD-TTD coarse grid J = 
  % PFD random J = 1.4863

  %          ╭──────────────────────────────────────────────────────────╮
    %          │                          plots                           │
    %          ╰──────────────────────────────────────────────────────────╯
  t = (0:F);
  tilefigure('Name','input')
  plot(t,u)
  grid on
  xlabel('Time step $k$','Interpreter','latex')
  ylabel('Input $\mathbf{u}_{k}$','Interpreter','latex')

  tilefigure('Name','statex')
  plot(t,x)
  grid on
  xlabel('Time step $k$','Interpreter','latex')
  ylabel('State $\mathbf{x}_{k}$','Interpreter','latex')

  tilefigure('Name','statemudecision')
  plot(t,mu,'o')
  grid on
  hold on
  plot(t,d,'x')
  xlabel('Time step $k$','Interpreter','latex')
  ylabel('State $\mu_{k}$, decision $d_{k}$','Interpreter','latex')
  legend({'model index $\mu_{k}$','decision $d_{k}$'},'Interpreter','latex')

  tilefigure
  plot(t,y)
  grid on
  xlabel('Time step $t_{k}$ [s]','Interpreter','latex')
  ylabel('Output $\mathbf{y}(t_{k})$','Interpreter','latex')



  for k = 1:Fp1
    Pxxf111(k) = estimate(k).Pxxfseq(1,1,1);
    Pxxf121(k) = estimate(k).Pxxfseq(1,2,1);
    Pxxf221(k) = estimate(k).Pxxfseq(2,2,1);
    Pxxf112(k) = estimate(k).Pxxfseq(1,1,2);
    Pxxf122(k) = estimate(k).Pxxfseq(1,2,2);
    Pxxf222(k) = estimate(k).Pxxfseq(2,2,2);
    pmuf(:,k) = estimate(k).pmufseq;
  end

  tcl = tiledlayout(3,2);
  nexttile
  plot(t,Pxxf111)
  grid on
  nexttile
  plot(t,Pxxf121)
  grid on
  nexttile
  plot(t,Pxxf221)
  grid on
  nexttile
  plot(t,Pxxf112)
  grid on
  nexttile
  plot(t,Pxxf122)
  grid on
  nexttile
  plot(t,Pxxf222)
  grid on
  title(tcl,'Variances')

  tilefigure
  plot(t,pmuf)
  grid on
  xlabel('Time step $t_{k}$ [s]','Interpreter','latex')
  ylabel('Output $\mathbf{y}(t_{k})$','Interpreter','latex')

  tilefigure
  plot(t(1:end-1),diff(x,[],2))
  ylabel('State differences')

end

tilefigure([nan nan])

