clear; close all; clc;
par_default;
    
%addpath(genpath('/home/gary/Documents/Simulations/'))
%% SET MODEL PARAMETERS

runID   = 'layer';  % run identifier tag
outdir  = '../out'; % output directory 
nout    = 10;       % print output every 'nout' steps
lvplt   = 1;        % plot figures live (1) or in background (0)     
svout   = 0;        % save figures and data to file (1)
bnchm   = 0;        % run benchmark

% set domain parameters
Nz      = 200;      % num. grid size
Nx      = Nz;
D       = 2e3;      % phys. domain depth [m]

% set physical parameters
mu      = 1e-4;     % pore fluid viscosity (water) [Pa s]
k0      = 1e-10;    % background permeability [m2]
n       = 3;        % permeability powerlaw [1]
rhol0   = 1000;     % fluid density [kg/m3]
grav    = 9.81;     % gravity [m/s2]
kC      = 1e-8;     % chemical diffusivity [m2/s]  
kT      = 1e-6;     % thermal diffusivity [m2/s]
kV      = 1e-7;     % vapour bubble diffusivity [m2/s]
aT      = 1e-4;     % thermal expansivity [1/K]
aC      = -1.1;     % chemical expansivity [1/wt]
aV      = 0.90;     % vapour-liquid density contrast [1/wt]
LH      = 2.25e6/4200; % temperature jump from latent heat of vapourisation [C]

% set initial condition parameters
finit   = 'layer';  % initial condition: 'linear' or 'layer'
f0      = 0.10;     % top/background initial porosity [vol]
f1      = f0;       % base porosity [vol]  
df      = 0.001;    % perturbation amplitude [vol]

Tinit   = 'layer';  % initial condition: 'linear' or 'layer'
Ttop    = 10;       % top boundary temperature
Tbot    = 100;      % base boundary temperature
T0      = Ttop;     % top/background initial temperature [C]
T1      = Tbot;     % base initial temperature [C]
dT      = -0.1;     % perturbation amplitude [C]

Cinit   = 'layer';  % initial condition: 'linear' or 'layer'
Ctop    = 0.010;    % top boundary concentration [wt]
Cbot    = 0.001;    % base boundary concentration [wt]
C0      = Ctop;     % top/background concentration  [wt]
C1      = Cbot;     % base concentration [wt]
dC      = 1e-4;     % perturbation amplitude [wt]

zlay    = 0.5;      % relative depth of layer boundary
wlay    = 0.02;     % relative width of layer boundary

smth    = 10;       % smoothness of initial fields

% set boundary conditions
BC_T    = {'closed','periodic'};
BC_C    = {'closed','periodic'};
BC_V    = {'closed','periodic'};
BC_VP   = {'closed','periodic'};

% set model timing parameters
tend    = 1e11;      % model stopping time [s]
Nt      = 1e4;       % max number of time step
dt      = 1e9;       % initial time step

% set numerical solver parameters
CFL     = 0.75;      % Courant-Friedrich-Lewy number to limit time step size
ADVN    = 'weno5';   % advection scheme
nup     = 100;       % update TC-solution and check residuals every nup iter
tol     = 1e-8;      % residual tolerance for iterative solver
maxit   = 1e4;       % maximum number of iterations
alpha   = 0.99;      % step size for iterative solver
beta    = 0.97;      % damping parameter for iterative solver

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************
