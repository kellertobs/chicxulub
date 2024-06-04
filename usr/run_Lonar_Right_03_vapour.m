clear; close all; clc;

par_Lonar_Right_03;

% SET MODEL PARAMETERS
runID   = 'Lonar_Right_03_vapour'; % run identifier tag
nout    = 20;       % print output every 'nop' steps
svout   = 1;        % save figures and data to file (1)

Nz      = 200;
Nx      = 200;

indir   = '../img_inputs/Lonar/Lonar_Right_03_200x200/'; % input directory for arrays


%% Make background temperature and porosity from an array (make array using ImgPrep_Temperature.m and ImgPrep_Porosity.m)
TArr   = load([indir 'Lonar_Right_01_200x200_TArray.mat']);
TArray = TArr.T_array2/2.25;

addpath ../src
for i = 1:300
    TArray = TArray + diffus(TArray,ones(size(TArray))/8,1,[1,2],BC_VP);
end

% fArr = load([indir 'Lonar_Right_01_200x200_fArray.mat']);
% fArray = fArr.f_array2;


%% Make lithologic units from binary images (2D arrays) prepared using ImgPrep_LithUnits.m
cbrc = imread(indir + "Lonar_Right_03_200x200_Lith_6.png"); % Binary image showing location of Coarse Breccia
sed  = imread(indir + "Lonar_Right_03_200x200_Lith_5.png"); % Binary image showing location of Sediments
bslt_pdr = imread(indir + "Lonar_Right_03_200x200_Lith_3.png"); % Binary image showing location of Powdered Basalt
mbrc = imread(indir + "Lonar_Right_03_200x200_Lith_8.png"); % Binary image showing location of Micro Breccia

bslt = imread(indir + "Lonar_Right_03_200x200_Lith_2.png"); % Binary image showing location of Basalt
bsmt = imread(indir + "Lonar_Right_03_200x200_Lith_1.png"); % Binary image showing location of Basement
wat  = imread(indir + "Lonar_Right_03_200x200_Lith_4.png"); % Binary image showing location of Water
air  = imread(indir + "Lonar_Right_03_200x200_Lith_7.png"); % Binary image showing location of Air

wat_evolve = false;              % evolve water as well-mixed reservoir; else keep T,C constant
tau_eqlb   = 1*3600*24*365.25;  % water-air thermal equilibration time 

wat(1,:)   = 0;                 % input image has wrong values along topmost row, please fix
air(1,:)   = 1;                 % input image has wrong values along topmost row, please fix

% [Xo,Zo]  = meshgrid(linspace(1,Nx,Nz),1:Nz);
% [Xi,Zi]  = meshgrid(1:Nx,1:Nz);
% cbrc     = logical(interp2(Xo,Zo,double(cbrc    ),Xi,Zi));
% sed      = logical(interp2(Xo,Zo,double(sed     ),Xi,Zi));
% bslt_pdr = logical(interp2(Xo,Zo,double(bslt_pdr),Xi,Zi));
% mbrc     = logical(interp2(Xo,Zo,double(mbrc    ),Xi,Zi));
% bslt     = logical(interp2(Xo,Zo,double(bslt    ),Xi,Zi));
% bsmt     = logical(interp2(Xo,Zo,double(bsmt    ),Xi,Zi));
% wat      = logical(interp2(Xo,Zo,double(wat     ),Xi,Zi));
% air      = logical(interp2(Xo,Zo,double(air     ),Xi,Zi));
% TArray   = interp2(Xo,Zo,TArray,Xi,Zi);

indstruct = cat(3,   air,   wat,   sed,   cbrc,   bslt_pdr,  mbrc,   bslt,   bsmt);   % Name and order of structures
fstruct   =       [f_air, f_wat, f_sed,  f_plb,      f_mlb, f_imr, f_bslt, f_bsmt];   % porosity of structures (nan = do not set)
Tstruct   =   nan*[T_air, T_wat, T_sed,  T_plb,      T_mlb, T_imr, T_bslt, T_bsmt];   % temperature of structures (nan = do not set)
Cstruct   =       [C_air, C_wat, C_sed,  C_plb,      C_mlb, C_imr, C_bslt, C_bsmt];   % salinity of structures (nan = do not set)

tol     = 1e-7;      % residual tolerance for iterative solver
alpha   = 1.25;      % step size for iterative solver
beta    = 0.99;      % damping parameter for iterative solver
maxit   = 5e3;       % maximum number of iterations

%*****  RUN CHICXULUB MODEL  *************************************************
run('../src/main')
%**************************************************************************
