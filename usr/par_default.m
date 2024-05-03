%% SET MODEL PARAMETERS

runID   = 'impact_default'; % run identifier tag
outdir  = '../out'; % output directory 
indir   = '../img_inputs/Kardla/Kardla_Right_200x200/'; % input directory for arrays
nout    = 50;       % print output every 'nop' steps
lvplt   = 1;        % plot figures live (1) or in background (0)     
svout   = 1;        % save figures and data to file (1)
bnchm   = 0;        % run benchmark

% set domain parameters
N       = 200;      % num. grid size
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


%% Make background temperature and porosity from an array (make array using ImgPrep_Temperature.m and ImgPrep_Porosity.m)
% TArr   = load([indir 'Kardla_Right_200x200_TArray.mat']);
% TArray = TArr.T_array2;
% 
% fArr = load([indir 'Kardla_Right_200x200_fArray.mat'])
% fArray = fArr.f_array2;
% 
% 
% %% Make lithologic units from binary images (2D arrays) prepared using ImgPrep_LithUnits.m
% svt = imread(indir + "Kardla_Right_200x200_Lith_4.png"); % Binary image showing location of Suevite
% plb = imread(indir + "Kardla_Right_200x200_Lith_3.png"); % Binary image showing location of Polymict Lithic Breccia
% sed = imread(indir + "Kardla_Right_200x200_Lith_5.png"); % Binary image showing location of Sediments
% wat = imread(indir + "Kardla_Right_200x200_Lith_2.png"); % Binary image showing location of Water
% 
% % Lithology parameters
% f_svt = 0.3;      % Suevite porosity
% T_svt = 600;      % Suevite Temperature
% C_svt = 0.01;     % Suevite salinity
% 
% f_plb = 0.2;      % Polymict Lithic Breccia porosity
% T_plb = 400;      % Polymict Lithic Breccia Temperature
% C_plb = 0.01;     % Polymict Lithic Breccia salinity
% 
% f_mlb = 0.1;      % Monomict Lithic Breccia porosity
% T_mlb = 400;      % Monomict Lithic Breccia Temperature
% C_mlb = 0.01;     % Monomict Lithic Breccia salinity
% 
% f_imr = 0.05;     % Impact Melt Rock porosity
% T_imr = 900;      % Impact Melt Rock Temperature
% C_imr = 0.01;     % Impact Melt Rock salinity
% 
% f_sed =  0.2;     % Sediment porosity (volume fraction) 
% T_sed =  100;     % Sediment Temperature
% C_sed =  0.01;     % Sediment salinity (wt fraction - 0.2 = 20% of mass is salinity_ sea water is 3.5% (0.035)
% 
% f_wat =  0.5;     % Water porosity
% T_wat =  100;     % Water Temperature
% C_wat =  0.035;   % Water salinity (wt fraction - 0.2 = 20% of mass is salinity_ sea water is 3.5% (0.035)
% 
% indstruct = cat(3,   wat,   plb,   svt,   sed);
% % fstruct   =       [f_wat, f_plb, f_svt, f_sed];   % porosity of structures (nan = do not set)
% fstruct   =       [  nan,   nan,   nan,   nan];   % porosity of structures (nan = do not set)
% % Tstruct   =       [T_wat, nan, T_svt, nan];      % temperature of structures (nan = do not set)
% Tstruct   =       [  nan,   nan,   nan,   nan];
% Cstruct   =       [C_wat, C_plb, C_svt, C_sed];   % salinity of structures (nan = do not set)
% watind = 1;
% Twat   = 10;

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
nup     = 100;       % update TC-solution and check residuals every nup iter
tol     = 1e-9;      % residual tolerance for iterative solver
maxit   = 1e4;       % maximum number of iterations
alpha   = 0.99;      % step size for iterative solver
beta    = 0.97;      % damping parameter for iterative solver
