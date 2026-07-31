clear; close all; clc;
par_default;

%addpath(genpath('/home/gary/Documents/Simulations/'))
%% SET MODEL PARAMETERS

runID   = 'Ra1e4'; % run identifier tag
outdir  = '../out'; % output directory 
nout    = 10;       % print output every 'nout' steps
lvplt   = 1;        % plot figures live (1) or in background (0)     
svout   = 1;        % save figures and data to file (1)
bnchm   = 0;        % run benchmark

% set domain parameters
Nz      = 200;      % num. grid size
Nx      = Nz;
D       = 1e3;      % phys. domain depth [m]

% set initial condition parameters
finit   = 'linear'; % initial condition: 'linear' or 'layer'
f0      = 0.20;     % top/background initial porosity [vol]
f1      = f0;       % base porosity [vol]  
df      = 0.0;      % perturbation amplitude [vol]

Tinit   = 'linear'; % initial condition: 'linear' or 'layer'
Ttop    = 10;       % top boundary temperature
Tbot    = 72.5;     % base boundary temperature
T0      = (Ttop+Tbot)/2; % top/background initial temperature [C]
T1      = (Ttop+Tbot)/2; % base initial temperature [C]
dT      =-T0/50;    % perturbation amplitude [C]

Cinit   = 'linear'; % initial condition: 'linear' or 'layer'
Ctop    = 0.001;    % top boundary concentration [wt]
Cbot    = 0.001;    % base boundary concentration [wt]
C0      = Ctop;     % top/background concentration  [wt]
C1      = Cbot;     % base concentration [wt]
dC      = 0;        % perturbation amplitude [wt]

aT      = 1e-4;     % thermal expansivity [1/K]
k0      = 2e-8;     % background permeability [m2]
grav    = 10;       % gravity [m/s2]
smth    = 5;        % smoothness of initial fields
bnd_w   = 25;       % initial boundary layer thickness

BDT     = 1e16;

% set boundary conditions
BC_T    = {[Ttop,Tbot],'periodic'};
BC_C    = {[Ctop,Cbot],'periodic'};
BC_V    = {'closed','periodic'};
BC_VP   = {'closed','periodic'};

% set model timing parameters
tend    = 1e11;      % model stopping time [s]
Nt      = 1e4;       % max number of time step
dt      = 1e9;       % initial time step

solver  = 'direct';  % select solver algorithm ('direct' or 'pseudotransient')
nup     = 1;         % update TC-solution and check residuals every nup iter
tol     = 1e-5;      % residual tolerance for iterative solver
maxit   = 10;        % maximum number of iterations

RaT     = (Tbot-Ttop)*aT*rhol0*grav*D*k0*f0^n/mul/kT


%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************
