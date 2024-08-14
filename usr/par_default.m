%% SET MODEL PARAMETERS

runID   = 'impact_default'; % run identifier tag
outdir  = '../out'; % output directory 
indir   = '../img_inputs/Kardla/Kardla_Right_200x200/'; % input directory for arrays
nout    = 50;       % print output every 'nop' steps
lvplt   = 1;        % plot figures live (1) or in background (0)     
svout   = 1;        % save figures and data to file (1)
bnchm   = 0;        % run benchmark

% set domain parameters
Nz      = 200;      % num. grid size
Nx      = 200;      % num. grid size
D       = 1e3;      % phys. domain depth [m]

% set physical parameters
mu      = 1e-4;     % pore fluid viscosity (water) [Pa s]
k0      = 1e-10;    % background permeability [m2]
n       = 3;        % permeability powerlaw [1]
rhol0   = 1000;     % fluid density [kg/m3]
grav    = 9.81;     % gravity [m/s2]
kC      = 1e-8;     % chemical diffusivity [m2/s]  
kT      = 1e-6;     % thermal diffusivity [m2/s]
kV      = 1e-10;    % vapour bubble diffusivity [m2/s]
aT      = 4e-4;     % thermal expansivity [1/K]
aC      = -0.7;     % chemical expansivity [1/wt]
aV      = 0.90;     % vapour-liquid density contrast [1/wt]
LH      = 2.25e6/4.2e3; % temperature jump from latent heat of vapourisation [C]

% set initial condition parameters
finit   = 'array'; % initial condition: 'linear' or 'layer' or 'array'
f0      = 0.10;     % top/background initial porosity [vol]
f1      = 0.01;     % base porosity [vol]  
df      = 0.001;    % perturbation amplitude [vol]

Tinit   = 'array'; % initial condition: 'linear' or 'layer' or 'array'
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

T_air      = 10;
T_wat      = 10;
C_wat      = 0.035;
wat_evolve = false;             % evolve water as well-mixed reservoir; else keep T,C constant
tau_eqlb   = 5*3600*24*365.25;  % water-air thermal equilibration time 
unit       = nan;

smth    = 10; % smoothness of initial fields

% set boundary conditions
BC_T    = {[Ttop,Tbot],'closed'};
BC_C    = {[Ctop,Cbot],'closed'};
BC_V    = {'closed','closed'};
BC_VP   = {'open','closed'};

% set model timing parameters
tend    = 1e11;      % model stopping time [s]
Nt      = 1e5;       % max number of time step
dt      = 1e3;       % initial time step

% set numerical solver parameters
CFL     = 1.0;       % Courant-Friedrich-Lewy number to limit time step size
ADVN    = 'weno5';   % advection scheme
nup     = 50;        % update TC-solution and check residuals every nup iter
tol     = 1e-8;      % residual tolerance for iterative solver
maxit   = 5e3;       % maximum number of iterations
alpha   = 0.99;      % step size for P-iterations
beta    = 0.97;      % damping parameter for P-iterations
gamma   = 0.25;      % step size for TC-iterations
delta   = 0.25;      % damping parameter for TC-iterations
