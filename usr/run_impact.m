clear; close all; clc;
par_default;

%% SET MODEL PARAMETERS

runID   = 'test';   % run identifier tag
outdir  = '../out'; % output directory 
indir   = '../img_inputs/Kardla/Kardla_Right_200x200/'; % input directory for arrays
nout    = 1;        % print output every 'nop' steps
lvplt   = 1;        % plot figures live (1) or in background (0)     
svout   = 0;        % save figures and data to file (1)
bnchm   = 0;        % run benchmark

% set domain parameters
Nz      = 200;      % num. grid size
Nx      = 300;
D       = 1e3;      % phys. domain depth [m]

% set physical parameters
mu      = 1e-4;     % pore fluid viscosity (water) [Pa s]
k0      = 1e-11;    % background permeability [m2]
n       = 3;        % permeability powerlaw [1]
rhol0   = 1000;     % fluid density [kg/m3]
grav    = 9.81;     % gravity [m/s2]
kC      = 1e-8;     % chemical diffusivity [m2/s]  
kT      = 1e-6;     % thermal diffusivity [m2/s]
aT      = 1e-4;     % thermal expansivity [1/K]
gC      = 1.1;      % chemical expansivity [1/wt]

% set initial condition parameters
finit   = 'linear'; % initial condition: 'linear' or 'layer' or 'array'
f0      = 0.10;     % top/background initial porosity [vol]
f1      = 0.01;     % base porosity [vol]  
df      = 0.001;    % perturbation amplitude [vol]

Tinit   = 'linear'; % initial condition: 'linear' or 'layer' or 'array'
Ttop    = 0;       % top boundary temperature
Tbot    = 40;       % base boundary temperature
T0      = 400;      % top/background initial temperature [C]
T1      = 40;       % base initial temperature [C]
dT      = -5;       % perturbation amplitude [C]

Cinit   = 'linear'; % initial condition: 'linear' or 'layer'
Ctop    = 0.035;    % top boundary concentration [wt]
Cbot    = 0.001;    % base boundary concentration [wt]
C0      = 0.010;    % top/background concentration  [wt]
C1      = 0.001;    % base concentration [wt]
dC      = 5.e-4;    % perturbation amplitude [wt]

zlay    = 0.5;      % relative depth of layer boundary
wlay    = 0.02;     % relative width of layer boundary

T_air   = 10;
T_wat   = 10;

smth    = 10; % smoothness of initial fields

% set boundary conditions
BC_T    = {[Ttop,Tbot],'closed'};
BC_C    = {[Ctop,Cbot],'closed'};
BC_VP   = {'open','closed'};

% set model timing parameters
tend    = 1e11;      % model stopping time [s]
Nt      = 1e4;       % max number of time step


% set numerical solver parameters
CFL     = 0.75;      % Courant-Friedrich-Lewy number to limit time step size
ADVN    = 'weno5';   % advection scheme
nup     = 50;       % update TC-solution and check residuals every nup iter
tol     = 1e-9;      % residual tolerance for iterative solver
maxit   = 1e4;       % maximum number of iterations
alpha   = 0.99;      % step size for iterative solver
beta    = 0.95;      % damping parameter for iterative solver

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************
