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
Nx      = Nz*1.5;
D       = 3e3;      % phys. domain depth [m]

% set permeability
k0      = 2e-10;    % background permeability [m2]
kD      = 1e-17;    % hot/ductile permeability [m2]
BDT     = 500;      % brittle/ductile transition T [C]
ld      = 1;        % dispersivity length [m2]

% set initial condition parameters
finit   = 'linear';  % initial condition: 'linear' or 'layer'
f0      = 0.25;     % top/background initial porosity [vol]
f1      = 0.05;     % base porosity [vol]  
df      = 0.001;    % perturbation amplitude [vol]

Tinit   = 'linear';  % initial condition: 'linear' or 'layer'
Ttop    = -25;       % top boundary temperature
Tbot    = 400;      % base boundary temperature
T0      = Ttop;     % top/background initial temperature [C]
T1      = Tbot;     % base initial temperature [C]
dT      =-Tbot/100; % perturbation amplitude [C]

Cinit   = 'linear';  % initial condition: 'linear' or 'layer'
Ctop    = 0.001;     % top boundary concentration [wt]
Cbot    = 0.001;    % base boundary concentration [wt]
C0      = Ctop;     % top/background concentration  [wt]
C1      = Cbot;     % base concentration [wt]
dC      = Cbot/100; % perturbation amplitude [wt]

zlay    = 0.75;     % relative depth of layer boundary
wlay    = 0.01;    % relative width of layer boundary
bnd_w   = 0;

smth    = 3;        % smoothness of initial fields

unit    = zeros(Nz,Nx);
unit(end-2*Nz/10:end-Nz/10+1,(Nx/10+1:end-Nx/10)) = 1;
funit   = 0.05;
Cunit   = 0.05;
Tunit   = 850;

T_air   = 10;
T_wat   = 10;
C_wat   = 0.001;

air     = zeros(Nz,Nx);
air(1:Nz/20,:) = 1;
wat     = zeros(Nz,Nx);
wat(Nz/20+1:2*Nz/20,:) = 1;
wat_evolve = 1;

% set boundary conditions
BC_T    = {'closed','periodic'};
BC_C    = {'closed','periodic'};
BC_V    = {'closed','periodic'};
BC_VP   = {'open','periodic'};

% set model timing parameters
tend    = 1e13;       % model stopping time [s]
Nt      = 1e4;        % max number of time step
dt      = 1e9;        % initial time step

% set numerical tuning parameters
CFL     = 0.75;
solver  = 'direct';  % select solver algorithm ('direct' or 'pseudotransient')
nup     = 1;         % update TC-solution and check residuals every nup iter
tol     = 1e-5;      % residual tolerance for iterative solver
maxit   = 10;        % maximum number of iterations
% solver  = 'pseudotransient'; % select solver algorithm ('direct' or 'pseudotransient')
% nup     = 10;         % update TC-solution and check residuals every nup iter
% tol     = 1e-6;        % residual tolerance for iterative solver
% maxit   = 1e3;         % maximum number of iterations

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************
