close all
clear all
clc

% Add path to functions needed in this example
addpath('../dttd')


scenario = 'second-order-oscillating';
valueFunctionFileName = 'bellman11D_finegrid_1times';
% valueFunctionFileName = 'bellman11D_finegrid';
simulationResultsFileName = 'simulationResults11D_finegrid_1times';
% simulationResultsFileName = 'simulationResults11D_finegrid';
FineGrid = true;
ProbFirst = false; % If probability is the first element

% trainBellman = 'fresh';
%trainBellman = 'continue';
% trainBellman = 'no';
trainBellman = 'fresh_11D_example';
% trainBellman = 'NN';

% simulationResultsFileName = 'simulationResults_passive';

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
  model.M(1).A = ssdata(c2d(ss(A,B,C,D),Ts));

  % Transition probability matrix
  model.P = [0.98 0.01
  0.02 0.99];

  % Detection cost function - not used
  Ld = [0 1
  1 0];

  % Discount factor
  eta = 0.9;

  % The predictive mean of the state
  meanx0p = [0
  0];

  % The predictive covariance of the state
  % Pxx0p = eye(2)*0.1; %
  Pxx0p = diag([7e-5 2e-4]);

  % Predictive probabilities of models
  pmu0p = [1
  0];

 %          ╭──────────────────────────────────────────────────────────╮
%          │     Define quantization levels for intermediate grid     │
%          ╰──────────────────────────────────────────────────────────╯

  if FineGrid
    nLvlMean = 61;
    nLvlvar = 10;
    nLvlProb = 1e3;

    xGridLvl1 = linspace(-1,1,nLvlMean); %5
    xGridLvl2 = linspace(-2,2,nLvlMean); %21
    xGridLvl3 = linspace(-0.5,0.5,nLvlMean); %11
    xGridLvl4 = linspace(-0.5,0.5,nLvlMean); %11
    ellGridLvl11 = sqrt(linspace(0.4e-5,0.6e-5,nLvlvar)); %4
    ellGridLvl22 = sqrt(linspace(2e-5,3e-5,nLvlvar)); %4
    anglesGridLvl = linspace(0,pi,10);
    % anglesGridLvl = [pi/2-0.1 pi/2 pi/2+0.1];
    probGridLvl = linspace(0,1,nLvlProb);
  else
    xGridLvl1 = -1:0.5:1;
    xGridLvl2 = -2:0.2:2;
    xGridLvl3 = -0.5:0.1:0.5;
    xGridLvl4 = -0.5:0.1:0.5;
    ellGridLvl11 = sqrt(linspace(0.4e-5,0.6e-5,4));
    ellGridLvl22 = sqrt(linspace(2e-5,3e-5,4));
    anglesGridLvl = [pi/2-0.1 pi/2 pi/2+0.1];
    probGridLvl = 0:0.01:1;
  end

  if ProbFirst
    thetaGridLvl = {probGridLvl,...
    xGridLvl1,xGridLvl2,xGridLvl3,xGridLvl4,...
    ellGridLvl11,...
    ellGridLvl22,anglesGridLvl,...
    ellGridLvl11,...
    ellGridLvl22,anglesGridLvl};
  else
    thetaGridLvl = {xGridLvl1,xGridLvl2,xGridLvl3,xGridLvl4,...
    ellGridLvl11,...
    ellGridLvl22,anglesGridLvl,...
    ellGridLvl11,...
    ellGridLvl22,anglesGridLvl,...
    probGridLvl};
  end
  thetaGridInnerEdge = generateinneredges(thetaGridLvl);

  uGridLvl = {-10:10:10};

%          ╭──────────────────────────────────────────────────────────╮
%          │                  Grid mapping functions                  │
%          ╰──────────────────────────────────────────────────────────╯
  gsubidx2theta = @(subidx) deaggregationidx(thetaGridLvl,subidx);
  % gtheta2xi = @(theta) [theta(1:4)
        %     sphericalcoor2covmatltvec(theta(5:7))
        %     sphericalcoor2covmatltvec(theta(8:10))
        %     theta(11)];
	if ProbFirst
		gtheta2xi = @(theta) [theta(2:3)
		theta(2:3)-theta(4:5)
		sphericalcoor2covmatltvec(theta(6:8))
		sphericalcoor2covmatltvec(theta(9:11))
		theta(1)];
	else
		gtheta2xi = @(theta) [theta(1:2)
		theta(1:2)-theta(3:4)
		sphericalcoor2covmatltvec(theta(5:7))
		sphericalcoor2covmatltvec(theta(8:10))
		theta(11)];
	end
%

%          ╭──────────────────────────────────────────────────────────╮
%          │          Inverses of the grid mapping functions          │
%          ╰──────────────────────────────────────────────────────────╯
        % gxi2theta = @(xi) [xi(1:4)
        %     covmatltvec2sphericalcoor(xi(5:7))
        %     covmatltvec2sphericalcoor(xi(8:10))
        %     xi(11)];
  if ProbFirst
    gxi2theta = @(xi) [xi(11)
    xi(1:2)
    xi(1:2)-xi(3:4)
    covmatltvec2sphericalcoor(xi(5:7))
    covmatltvec2sphericalcoor(xi(8:10))];
  else
    gxi2theta = @(xi) [xi(1:2)
    xi(1:2)-xi(3:4)
    covmatltvec2sphericalcoor(xi(5:7))
    covmatltvec2sphericalcoor(xi(8:10))
    xi(11)];
  end

otherwise
  error('Unknown scenario')
end

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
if ProbFirst
  gtheta2subidx_control = @(theta) aggregationidx([theta(2:11);theta(1)],thetaGridInnerEdge); % picks nearest neighbor
else
  gtheta2subidx_control = @(theta) aggregationidx(theta,thetaGridInnerEdge); % picks nearest neighbor
end
% gtheta2subidx = @(theta) aggregationidx(theta,thetaGridInnerEdge); % picks nearest neighbor
% gtheta2subidx = @(theta) aggregationidx_reg_grid(theta,thetaGridInnerEdge); % works only for regular grids
gtheta2subidx = @(theta) theta2realidx(theta,thetaGridInnerEdge); % real idx

% Define system wraped into grid mappings
fsux_h = @(subidx,u,y) gtheta2subidx(gxi2theta(phismmkf(gtheta2xi(gsubidx2theta(subidx)),u,y,model,vec2tril,tril2vec,reducedProbabilityVector)));

% Define predicitve distribution of future output
pxsugs_h = @(subidx,u) mcwynext(gtheta2xi(gsubidx2theta(subidx)),u,model,vec2tril,reducedProbabilityVector);

% Define cost functions
mL_h = @(subidx,u)  min(probGridLvl(subidx(end)),1-probGridLvl(subidx(end)));
%mV_h = @(subidx,u,V) mVutttd(subidx,u,V,fsux_h,pxsugs_h);
mV_h = @(subidx,u,V) mVutttd_lin_interp(subidx,u,V,fsux_h,pxsugs_h);

switch trainBellman
case 'fresh'
  V = vvittd(thetaGridLvl,uGridLvl,eta,mL_h,mV_h,[]);
  save('bellman','V')
case 'continue'
  load('bellman','V')
  V = vvittd(thetaGridLvl,uGridLvl,eta,mL_h,mV_h,V);
  save('bellman','V')
case 'fresh_11D_example'
  % [V,dVfrobnorm,dVfrobnorm_ci,dVfrobnorm_li] = vvittd_11D_example(thetaGridLvl,thetaGridInnerEdge);
        % [V,dVfrobnorm,dVfrobnorm_ci,dVfrobnorm_li] = vvittd_11D_example_ci(thetaGridLvl,thetaGridInnerEdge);
        % save('bellman31_ci','dVfrobnorm_li','dVfrobnorm_ci','thetaGridLvl')
  gridElements = prod(cellfun(@length,thetaGridLvl));
  if ProbFirst
    [V,dVfrobnorm1,dVfrobnorm2,dVfrobnorm3] = vvittd_11D_example_p1st(thetaGridLvl,thetaGridInnerEdge);
    save(valueFunctionFileName,'V','dVfrobnorm1','dVfrobnorm2','dVfrobnorm3','gridElements','thetaGridLvl')
  else
    [V,dVfrobnorm1,dVfrobnorm2,dVfrobnorm3] = vvittd_11D_example(thetaGridLvl,thetaGridInnerEdge);
    save(valueFunctionFileName,'V','dVfrobnorm1','dVfrobnorm2','dVfrobnorm3','gridElements','thetaGridLvl')
  end
case 'no'
  load(valueFunctionFileName)
case 'NN'
  % load('NN2_parameters_continuous');
  load('NN2_parameters_discrete_attempt2','problem','model');

  % problem = 'active_detec_estim_secondorder';
  % model = configModel(problem);
  % policy = configPolicy(model,problem);
  load('results_active_detec_2ndOrder_20240325T114121.mat','policy','PHI_star')
  x_dim = model.x_dim;
  u_dim = 1;
  neurons = policy.approx.neurons;
  theta = PHI_star.mu;
  W = reshape(theta(1:(x_dim+1)*neurons),(x_dim+1),neurons);
  C = reshape(theta((x_dim+1)*neurons+1:end),neurons,u_dim);
  % for continuous bounded policy
  % policy_NN = @(X) max(model.min_u,min(model.max_u,mlp_out(C,W,X,model,policy)));
  % for discrete policy
   % policy_NN = @(X) round(max(model.min_u,min(model.max_u,mlp_out(C,W,X,model,policy)))/model.max_u)*model.max_u;

  utmp = @(X) mlp_out(C,W,X,model,policy);
otherwise
  if ProbFirst
    load('./bellman_P1st','V')
  else
    % load('bellman','V')
  % load('./bellman_orig.mat','V')
    % load('bellman_finegrid_angle10_P1st')
    % load('./bellman_P1st','V')
    % load('./bellman_orig_e123.mat','V')
  % load('./bellman_50iter.mat','V')
    % load('./bellman_50iter_finegrid.mat','V')
    % load('./bellman_30iter_finegrid_10angles_e123.mat','V')
    % load('./bellman41.mat','V')
  % load('./bellman51.mat','V')
  % load('./bellman61.mat','V')
  % load('./bellman61_10_1000.mat','V')
  end
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

controller = @(xi) getvp(gtheta2subidx_control(gxi2theta(xi)),V,inputs,eta,mL_h,mV_h);
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
    for k=1:Fp1
      estimate(k).xpseq = zeros(nx,nModel);
      estimate(k).Pxxpseq = zeros(nx,nx,nModel);
      estimate(k).pmupseq = zeros(nModel,1);
      estimate(k).ypseq = zeros(ny,nModel);
      estimate(k).Pyypseq = zeros(ny,ny,nModel);
      estimate(k).xfseq = zeros(nx,nModel);
      estimate(k).Pxxfseq = zeros(nx,nx,nModel);
      estimate(k).pmufseq = zeros(nModel^2,1);
    end


    % Initialize estimator
    estimate(1).xpseq = meanx0p;
    estimate(1).Pxxpseq = Pxx0p;
    estimate(1).pmupseq = pmu0p;

    % generate initial state x and mu
    x(:,1) = meanx0p + chol(Pxx0p,'lower')*randn(nx,1);
    mu(1) = gendrnd(pmu0p);
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
      case {'no','fresh', 'continue', 'fresh_11D_example'}
        % form xi
        xi = sufstat2hyperstate(estimate(k).xfseq,estimate(k).Pxxfseq,estimate(k).pmufseq,tril2vec,reducedProbabilityVector);
        % Generate input
        [~,u(:,k)] = controller(xi);
        %%%%%%%%%%%%%%%%%%%%%%%% PFD
      case 'passive'
        u(:,k) = 10*randi([-1 1],1,1); % for passive FD
      % freq = 0.5;%[Hz]
      % Ts = 0.05; %same as for original system
      % u(:,k) = 10*cos(2*pi*freq*Ts*k);
      % u(:,k) = 0;
      %%%%%%%%%%%%%%%%%%%%%%%% NN AFD
      case 'NN'
        xi = sufstat2hyperstate(estimate(k).xfseq,estimate(k).Pxxfseq,estimate(k).pmufseq,model.tril2vec,model.reducedProbabilityVector);
        % u(:,k) = policy_NN(xi);
        u_ = utmp(xi);
        u_set = [-10,0,10];
        [~,jj] = min(dist(u_,u_set));
        u(:,k) = u_set(jj);
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
    J(iMC) = sum(eta.^t.*sum(d~=mu,1));
    fprintf('%u/%u MC simulations done.\n',iMC,nMC);
  end
  % save('simulationResults_finegrid_P1st','J');
  % save('simulationResults_NN_attempt2','J');
  save(simulationResultsFileName,'J');
% FineGrid P 1st 0.8406
% FineGrid P 1st angles10 1.3031
% Coarse grid P 1st angles3 0.8584

% FineGrid  0.8597
% Coarse grid 0.8683
% PFD 1.0650

% Adaptive computations (5 recomputations)
% Coarse grid 0.8264
% Fine grid 0.7574

% Quadratic interpolation Adaptive computations (1 recomputations)
% Coarse grid 0.7221
% Fine grid   0.7197
% Quadratic interpolation Adaptive computations (5 recomputations)
% Coarse grid 0.7404
% Fine grid 0.7678

% NN = 0.7664 - continuous u <-10,10>
% NN = 0.8482 - discrete u {-10,0,10}

% Coarse grid P1st 0.8253
% Fine grid P1st 0.9205
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

