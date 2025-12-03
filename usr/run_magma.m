clear; close all; clc;
par_default;
    
%% SET MODEL PARAMETERS

runID   = 'magma';  % run identifier tag
outdir  = '../out'; % output directory 
nout    = 20;       % print output every 'nout' steps
lvplt   = 1;        % plot figures live (1) or in background (0)     
svout   = 0;        % save figures and data to file (1)
bnchm   = 0;        % run benchmark

% set domain parameters
Nz      = 300;      % num. grid size
Nx      = Nz;
D       = 1e3;      % phys. domain depth [m]

% set permeability
k0      = 1e-12;    % background permeability [m2]
kD      = 1e-18;

% set initial condition parameters
finit   = 'layer';  % initial condition: 'linear' or 'layer'
f0      = 0.20;     % top/background initial porosity [vol]
f1      = 0.01;     % base porosity [vol]  
df      = 0.001;    % perturbation amplitude [vol]

Tinit   = 'layer';  % initial condition: 'linear' or 'layer'
Ttop    = 50;       % top boundary temperature
Tbot    = 1000;     % base boundary temperature
T0      = Ttop;     % top/background initial temperature [C]
T1      = Tbot;     % base initial temperature [C]
dT      =-Tbot/100; % perturbation amplitude [C]

Cinit   = 'layer';  % initial condition: 'linear' or 'layer'
Ctop    = 0.001;    % top boundary concentration [wt]
Cbot    = 0.020;    % base boundary concentration [wt]
C0      = Ctop;     % top/background concentration  [wt]
C1      = Cbot;     % base concentration [wt]
dC      = Cbot/100; % perturbation amplitude [wt]

zlay    = 0.75;     % relative depth of layer boundary
wlay    = 0.01;     % relative width of layer boundary

smth    = 2;        % smoothness of initial fields


% set boundary conditions
BC_T    = {'closed','periodic'};
BC_C    = {'closed','periodic'};
BC_V    = {'closed','periodic'};
BC_VP   = {'closed','periodic'};

% set model timing parameters
tend    = 1e11;      % model stopping time [s]
Nt      = 1e4;       % max number of time step
dt      = 1e9;       % initial time step

% set numerical tuning parameters
nup     = 100;       % update TC-solution and check residuals every nup iter
tol     = 1e-7;      % residual tolerance for iterative solver
maxit   = 5e3;       % maximum number of iterations
alpha   = 0.99;      % step size for P-iterations
beta    = 0.95;      % damping parameter for P-iterations
gamma   = 0.10;      % step size for TC-iterations
delta   = 0.00;      % damping parameter for TC-iterations

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************
