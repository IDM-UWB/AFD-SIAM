function [model,problem] = four_2nd_order_oscillating_models(reducedProbabilityVector)

           % ╭───────────────────────────────────────────────────────────╮
           % │                       DEFINE MODELS                       │
           % ╰───────────────────────────────────────────────────────────╯

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


  % The dimension of state and input
  [model.nx,model.nu] = size(model.M(1).B);

  % The dimension of measurement
  model.ny = size(model.M(1).C,1);

  % The dimension of state noise
  model.nw = size(model.M(1).G,2);

  % The dimension of measurement noise
  model.nv = size(model.M(1).H,2);

  % The number of models
  model.nModel = length(model.M);

  % Auxiliary indexing matrices
  model.tril2vec = idxtril2vec([model.nx,model.nx,model.nModel]);
  model.vec2tril = idxvec2tril([model.nx,model.nx,model.nModel]);

           % ╭───────────────────────────────────────────────────────────╮
           % │                      DEFINE PROBLEM                       │
           % ╰───────────────────────────────────────────────────────────╯
  FineGrid = false; 
  % Detection cost function - not used
  problem.Ld = ones(model.nModel) - eye(model.nModel);

  % Discount factor
  problem.eta = 0.9;

  if reducedProbabilityVector
    problem.nxi = model.nModel*(model.nx+model.nx*(model.nx+1)/2+1)-1;
  else
    problem.nxi = model.nModel*(model.nx+model.nx*(model.nx+1)/2+1);
  end
  %          ╭──────────────────────────────────────────────────────────╮
  %          │     Define quantization levels for intermediate grid     │
  %          ╰──────────────────────────────────────────────────────────╯

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
    % probGridLvl = linspace(1e-3,1-1e-3,100);
    probGridLvl = linspace(1e-3,1-1e-3,10);
  end

  problem.thetaGridLvl = {xGridLvl1,xGridLvl2,xGridLvl3,xGridLvl4,xGridLvl3,xGridLvl4,xGridLvl3,xGridLvl4,...
  ellGridLvl11,...
  ellGridLvl22,anglesGridLvl,...
  ellGridLvl11,...
  ellGridLvl22,anglesGridLvl,...
  ellGridLvl11,...
  ellGridLvl22,anglesGridLvl,...
  ellGridLvl11,...
  ellGridLvl22,anglesGridLvl,...
  probGridLvl,probGridLvl,probGridLvl};
  problem.thetaGridInnerEdge = generateinneredges(problem.thetaGridLvl);

  problem.sizV = cellfun(@length,problem.thetaGridLvl);
  problem.gridElements = prod(problem.sizV);

  problem.uGridLvl = {-10:10:10};

  %          ╭──────────────────────────────────────────────────────────╮
  %          │                  Grid mapping functions                  │
  %          ╰──────────────────────────────────────────────────────────╯
  problem.gsubidx2theta = @(subidx) deaggregationidx(problem.thetaGridLvl,subidx);
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
  problem.gtheta2subidx = @(theta) theta2realidx(theta,problem.thetaGridInnerEdge); % real idx

  % Define cost functions
  problem.mL_h = @(subidx,u)  min([probGridLvl(subidx(end-2:end)) 1-sum(probGridLvl(subidx(end-2:end)))]);
