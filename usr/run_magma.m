clear; close all; clc;
par_default;
    
%% SET MODEL PARAMETERS

runID   = 'magma';  % run identifier tag
outdir  = '../out'; % output directory 
nout    = 50;       % print output every 'nout' steps
lvplt   = 1;        % plot figures live (1) or in background (0)     
svout   = 1;        % save figures and data to file (1)
bnchm   = 0;        % run benchmark

% set domain parameters
Nz      = 200;      % num. grid size
Nx      = Nz;
D       = 1e3;      % phys. domain depth [m]

% set permeability
k0      = 1e-12;    % background permeability [m2]
kD      = 1e-18;    % hot/ductile permeability [m2]
BDT     = 500;      % brittle/ductile transition T [C]

% set initial condition parameters
finit   = 'linear';  % initial condition: 'linear' or 'layer'
f0      = 0.25;     % top/background initial porosity [vol]
f1      = 0.05;     % base porosity [vol]  
df      = 0.001;    % perturbation amplitude [vol]

Tinit   = 'linear';  % initial condition: 'linear' or 'layer'
Ttop    = 50;       % top boundary temperature
Tbot    = 400;      % base boundary temperature
T0      = Ttop;     % top/background initial temperature [C]
T1      = Tbot;     % base initial temperature [C]
dT      =-Tbot/100; % perturbation amplitude [C]

Cinit   = 'linear';  % initial condition: 'linear' or 'layer'
Ctop    = 0.001;    % top boundary concentration [wt]
Cbot    = 0.020;    % base boundary concentration [wt]
C0      = Ctop;     % top/background concentration  [wt]
C1      = Cbot;     % base concentration [wt]
dC      = Cbot/100; % perturbation amplitude [wt]

zlay    = 0.75;     % relative depth of layer boundary
wlay    = 0.01;    % relative width of layer boundary
bnd_w   = 0;

smth    = 2;        % smoothness of initial fields

unit    = zeros(Nz,Nx);
unit(end-30:end,:) = 1;
funit   = 0.04;
Cunit   = 0.03;
Tunit   = 800;

air     = zeros(Nz,Nx);
air(1:10,:) = 1;
wat     = zeros(Nz,Nx);
wat(11:20,:) = 1;
wat_evolve   = 1;
C_wat   = 0.0;

% set boundary conditions
BC_T    = {'closed','periodic'};
BC_C    = {'closed','periodic'};
BC_V    = {'closed','periodic'};
BC_VP   = {'closed','periodic'};

% set model timing parameters
tend    = 1e11;       % model stopping time [s]
Nt      = 1e4;        % max number of time step
dt      = 1e9;        % initial time step

% set numerical tuning parameters
nup     = 20;        % update TC-solution and check residuals every nup iter
tol     = 1e-8;       % residual tolerance for iterative solver
maxit   = 5e3;        % maximum number of iterations

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************
