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
D       = 1e3;      % phys. domain depth [m]

% set initial condition parameters
finit   = 'layer';  % initial condition: 'linear' or 'layer'
f0      = 0.20;     % top/background initial porosity [vol]
f1      = f0;       % base porosity [vol]  
df      = 0.0;      % perturbation amplitude [vol]

Tinit   = 'layer';  % initial condition: 'linear' or 'layer'
Ttop    = 10;       % top boundary temperature
Tbot    = 100;      % base boundary temperature
T0      = Ttop;     % top/background initial temperature [C]
T1      = Tbot;     % base initial temperature [C]
dT      =-Tbot/100; % perturbation amplitude [C]

Cinit   = 'layer';  % initial condition: 'linear' or 'layer'
Ctop    = 0.001;    % top boundary concentration [wt]
Cbot    = 0.020;    % base boundary concentration [wt]
C0      = Ctop;     % top/background concentration  [wt]
C1      = Cbot;     % base concentration [wt]
dC      = Cbot/100; % perturbation amplitude [wt]

zlay    = 0.5;      % relative depth of layer boundary
wlay    = 0.02;     % relative width of layer boundary

smth    = 5;        % smoothness of initial fields

% set boundary conditions
BC_T    = {'closed','periodic'};
BC_C    = {'closed','periodic'};
BC_V    = {'closed','periodic'};
BC_VP   = {'closed','periodic'};

% set model timing parameters
tend    = 1e11;      % model stopping time [s]
Nt      = 1e4;       % max number of time step
dt      = 1e9;       % initial time step


%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************
