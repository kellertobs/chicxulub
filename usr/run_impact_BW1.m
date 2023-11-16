clear; close all; clc;
    
%addpath(genpath('/home/gary/Documents/Simulations/'))
%% SET MODEL PARAMETERS

runID   = 'impact'; % run identifier tag
outdir  = '../out'; % output directory 
nout    = 50;       % print output every 'nop' steps
lvplt   = 1;        % plot figures live (1) or in background (0)     
svout   = 1;        % save figures and data to file (1)



% set domain parameters
N       = 200;      % num. grid size
D       = 2.5e3;      % phys. domain depth [m]
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
finit   = 'linear'; % initial condition: 'linear' or 'layer'
f0      = 0.20;     % top/background initial porosity [vol]
f1      = 0.005;    % base porosity [vol]  
df      = 0.001;    % perturbation amplitude [vol]

Tinit   = 'image'; % initial condition: 'linear' or 'layer'
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
whos

% % Trying to make background T from an image*********************
TempArray = imread("GreyScaleT_202x202_2.png");
TempArray = double(TempArray);
TempArray = TempArray.*3.33333;



% % LITHOLOGIC UNIT INPUTS***************************
svt = imread("Kardla_TopRight_202x202_Cluster 4.png"); % Binary image showing location of Suevite
plb = imread("Kardla_TopRight_202x202_Cluster 3.png"); % Binary image showing location of Polymict Lithic Breccia
mlb = imread("Kardla_TopRight_202x202_Cluster 2.png"); % Binary image showing location of Monomict Lithic Breccia
imr = imread("bwTest2.png"); % Binary image showing location of Impact Melt Rock
sed = imread("bwTest5.png"); % Binary image showing location of Sediments
wat = imread("Kardla_TopRightCluster 1_inv.png"); % Binary image showing location of Water


% Rocktype parameters
f_svt = 0.2;      % Suevite porosity
T_svt = 600;      % Suevite Temperature
C_svt = 0.01;     % Suevite salinity

f_plb = 0.1;      % Polymict Lithic Breccia porosity
T_plb = 400;      % Polymict Lithic Breccia Temperature
C_plb = 0.01;     % Polymict Lithic Breccia salinity

f_mlb = 0.1;      % Monomict Lithic Breccia porosity
T_mlb = 400;      % Monomict Lithic Breccia Temperature
C_mlb = 0.01;     % Monomict Lithic Breccia salinity

f_imr = 0.05;     % Impact Melt Rock porosity
T_imr = 900;      % Impact Melt Rock Temperature
C_imr = 0.01;      % Impact Melt Rock salinity

f_sed =  0.2;     % Sediment porosity (volume fraction) 
T_sed =  100;     % Sediment Temperature
C_sed = 0.01;     % Sediment salinity (wt fraction - 0.2 = 20% of mass is salinity_ sea water is 3.5% (0.035)

f_wat =  0.2;     % Water porosity
T_wat =  100;     % Water Temperature
C_wat =  0.035;     % Water salinity


% xstruct = [D/2,D/5,4*D/5];    % midpoint x-position of structures
% zstruct = [250,150,150];      % midpoint z-position of structures
% hstruct = [500,400,400];      % height of structures
% wstruct = [D,D/20,D/20];      % width of structures
% astruct = [0,30,-30];       % angle of structures to horizontal (counter-clockwise)

indstruct = cat(4,   wat,   plb,   svt,   sed);
fstruct   =       [f_wat, f_plb, f_svt, f_sed];   % porosity of structures (nan = do not set)
% Tstruct   =       [T_wat, T_plb, T_svt, T_sed];      % temperature of structures (nan = do not set)
Tstruct   =       [  nan,   nan,   nan,   nan];
Cstruct   =       [C_wat, C_plb, C_svt, C_sed];   % salinity of structures (nan = do not set)

smth    = 5*(N/100)^2; % smoothness of initial fields

% set boundary conditions
BC_T    = {[Ttop,Tbot],'closed'};
BC_C    = {[Ctop,Cbot],'closed'};
BC_VP   = {'open','closed'};

% set model timing parameters
tend    = 1e11;      % model stopping time [s]

% set numerical solver parameters
CFL     = 0.50;      % Courant-Friedrich-Lewy number to limit time step size
ADVN    = 'weno5';   % advection scheme
nup     = 200;       % update TC-solution and check residuals every nup iter
tol     = 1e-8;      % residual tolerance for iterative solver
maxit   = 1e4;       % maximum number of iterations
alpha   = 0.99;      % step size for iterative solver
beta    = 0.97;      % damping parameter for iterative solver

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main_comments.m')
%**************************************************************************
