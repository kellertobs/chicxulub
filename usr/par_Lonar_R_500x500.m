par_default;

%% SET MODEL PARAMETERS
%% Overall parameters (for run_impact and ImgPrep_ m-files)

runID   = 'Lonar_R'; % run identifier tag
outdir  = '../out'; % output directory 

% set domain parameters
Nz      = 500;      % num. grid size in z-direction
Nx      = 500;       % num. grid size in z-direction
D       = 1e3;      % phys. domain depth [m]


%% Set conditions for Img_Prep files
outdir_ImgPrep = '../img_inputs/Lonar/';

imgName_Lith   = "../img_inputs/Lonar/Lonar_Full_Lithology.png";
nUnits_Lith    = 8;

imgName_T      = "../img_inputs/Lonar/Lonar_Full_Temperature.png";  % Filename of the temperature distribution
nUnits_T       = 8;      % Need to specify number of Temperature units here

% imgName_f      = "../img_inputs/Lonar/Lonar_Right_Porosity.png"     % Filename of the image of the cross section
% nUnits_f       = 11;     % Need to specify number of Temperature units here; think they're subsequently clustered with 1 being highest area and 5 being lowest

% Crop function to select only part of the image defined by rectangle with [xmin ymin width height] REMEMBER: Origin is in top left for MatLab reasons
% These are scaling factors for splitting the image as fractions of width and height of the image. 
% Some basic options:
%     - Full image:         (0, 0, 1, 1)
%     - Top right of image: (0.5, 0, 0.5, 0.5)
%     - Top left of image:  (0, 0, 0.5, 0.5)

x_crp = 0.5; 
y_crp = 0;
w_crp = 0.5;
h_crp = 0.5;


%% Set conditions for run_impact.m

% set initial condition parameters
finit   = 'linear'; % initial condition: 'linear' or 'layer' or 'array'
f0      = 0.10;     % top/background initial porosity [vol]
f1      = 0.01;     % base porosity [vol]  
df      = 0.001;    % perturbation amplitude [vol]

Tinit   = 'linear'; % initial condition: 'linear' or 'layer' or 'array'
Ttop    = 0;        % top boundary temperature
Tbot    = 30;       % base boundary temperature
T0      = Ttop;     % top/background initial temperature [C]
T1      = Tbot;     % base initial temperature [C]
dT      = -1;       % perturbation amplitude [C]

Cinit   = 'linear'; % initial condition: 'linear' or 'layer'
Ctop    = 0.001;    % top boundary concentration [wt]
Cbot    = 0.01;     % base boundary concentration [wt]
C0      = Ctop;     % top/background concentration  [wt]
C1      = Cbot;     % base concentration [wt]
dC      = 5.e-4;    % perturbation amplitude [wt]


% Lithology parameters
f_svt = 0.3;      % Suevite porosity
T_svt = 600;      % Suevite Temperature
C_svt = 0.0;      % Suevite salinity

f_plb = 0.2;      % Polymict Lithic Breccia porosity
T_plb = 400;      % Polymict Lithic Breccia Temperature
C_plb = 0.0;      % Polymict Lithic Breccia salinity

f_mlb = 0.1;      % Monomict Lithic Breccia porosity
T_mlb = 400;      % Monomict Lithic Breccia Temperature
C_mlb = 0.0;      % Monomict Lithic Breccia salinity

f_imr = 0.05;     % Impact Melt Rock porosity
T_imr = 900;      % Impact Melt Rock Temperature
C_imr = 0.0;      % Impact Melt Rock salinity ***************************

f_sed =  0.25;    % Sediment porosity (volume fraction) 
T_sed =  100;     % Sediment Temperature
C_sed =  0.0;     % Sediment salinity (wt fraction)

f_wat =  1.0;     % Water porosity
T_wat =  10;      % Water Temperature (set to nan to evolve as mixed reservoir)
C_wat =  0.035;   % Water salinity (wt fraction, sea water is 3.5% (0.035), set to nan to evolve as mixed reservoir)

f_air =  0.0;     % Air porosity
T_air =  10;      % Air Temperature (will remain constant)
C_air =  0.0;     % Air salinity (wt fraction, set to 0, will remain constant)

f_bslt = 0.07;    % Basalt porosity************************************* 0.07 works, 0.05 doesn't
T_bslt = 100;     % Basalt Temperature 
C_bslt = 0.0;     % Basalt salinity

f_bsmt = 0.05;    % Basement porosity
T_bsmt = 200;     % Basement Temperature 
C_bsmt = 0.0;     % Basement salinity 


smth    = 5; % smoothness of initial fields

% set boundary conditions
BC_T    = {{'flux',[0,0.03]},'closed'};
BC_C    = {'closed','closed'};
BC_V    = {'closed','closed'};
BC_VP   = {'closed','closed'};
