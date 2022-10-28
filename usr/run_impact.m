clear; close all; clc;
    
%addpath(genpath('/home/gary/Documents/Simulations/'))
%% SET MODEL PARAMETERS

runID   = 'impact'; % run identifier tag
nop     = 50;       % print output every 'nop' steps
lvplt   = 1;        % plot figures live (1) or in background (0)     
svfig   = 1;        % save figures to file (1)

% set domain parameters
N       = 200;      % num. grid size
D       = 1e3;      % phys. domain depth [m]

% set physical parameters
mu      = 1e-3;     % pore fluid viscosity (water) [Pa s]
a       = 1e-3;     % grain size of matrix (sandstone) [m]
b       = 500;      % geom. factor for permeability [1]
n       = 3;        % permeability powerlaw [1]
rhol0   = 1000;     % fluid density [kg/m3]
grav    = 9.81;     % gravity [m/s2]
kC      = 1e-8;     % chemical diffusivity [m2/s]  
kT      = 1e-6;     % thermal diffusivity [m2/s]
aT      = 1e-4;     % thermal expansivity [1/K]
gC      = 1;        % chemical expansivity [1/wt]

% set initial condition parameters
finit   = 'linear'; % initial condition: 'linear' or 'layer'
f0      = 0.15;     % top/background initial porosity [vol]
f1      = 0.01;     % base porosity [vol]  
df      = 0.002;    % perturbation amplitude [vol]

Tinit   = 'linear'; % initial condition: 'linear' or 'layer'
Ttop    = 10;       % top boundary temperature
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

xstruct = [D/2];    % position of structures
zstruct = [0];      % onset depth of structures
dstruct = [500];    % final depth of structures
wstruct = [500];    % thickness of structures
astruct = [90];     % angle of structures to vertical
fstruct = [0.20];   % porosity of structures (nan = do not set)
Tstruct = [400];    % temperature of structures (nan = do not set)
Cstruct = [0.01];   % salinity of structures (nan = do not set)

smth    = 5*(N/100)^2; % smoothness of initial fields

% set boundary conditions
BC_T    = {[Ttop,Tbot],'periodic'};
BC_C    = {[Ctop,Cbot],'periodic'};
BC_VP   = {'open','periodic'};

% set model timing parameters
tend    = 1e11;      % model stopping time [s]

% set numerical solver parameters
CFL     = 0.50;      % Courant-Friedrich-Lewy number to limit time step size
ADVN    = 'weno5';   % advection scheme
nup     = 200;       % update TC-solution and check residuals every nup iter
tol     = 1e-8;      % residual tolerance for iterative solver
maxit   = 1e4;       % maximum number of iterations
alpha   = 0.99;      % step size for iterative solver
beta    = 0.98;      % damping parameter for iterative solver

% create output directory
if ~isfolder(['../out/',runID])
    mkdir(['../out/',runID]);
end

% run code
addpath ../src
main
