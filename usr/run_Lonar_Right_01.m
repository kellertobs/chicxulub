clear; close all; clc;

par_Lonar_Right_01;

% SET MODEL PARAMETERS
% runID   = 'Lonar'; % run identifier tag
nout    = 10;       % print output every 'nop' steps
svout   = 1;        % save figures and data to file (1)

indir   = '../img_inputs/Lonar/Lonar_Right_01_200x200/'; % input directory for arrays


%% Make background temperature and porosity from an array (make array using ImgPrep_Temperature.m and ImgPrep_Porosity.m)
TArr   = load([indir 'Lonar_Right_01_200x200_TArray.mat']);
TArray = TArr.T_array2;

% fArr = load([indir 'Lonar_Right_01_200x200_fArray.mat']);
% fArray = fArr.f_array2;


%% Make lithologic units from binary images (2D arrays) prepared using ImgPrep_LithUnits.m
bslt = imread(indir + "Lonar_Right_01_200x200_Lith_3.png"); % Binary image showing location of Basalt
bsmt = imread(indir + "Lonar_Right_01_200x200_Lith_1.png"); % Binary image showing location of Basement
wat = imread(indir + "Lonar_Right_01_200x200_Lith_4.png"); % Binary image showing location of Water (or Air)
cbrc = imread(indir + "Lonar_Right_01_200x200_Lith_6.png"); % Binary image showing location of Coarse Breccia
sed = imread(indir + "Lonar_Right_01_200x200_Lith_5.png"); % Binary image showing location of Sediments
bslt_pdr = imread(indir + "Lonar_Right_01_200x200_Lith_2.png"); % Binary image showing location of Powdered Basalt
mbrc = imread(indir + "Lonar_Right_01_200x200_Lith_7.png"); % Binary image showing location of Micro Breccia

indstruct = cat(3,   wat,   sed,   cbrc,   bslt_pdr,  mbrc,   bslt,   bsmt);   % Name and order of structures
fstruct   =       [f_wat, f_sed,  f_plb,      f_mlb, f_imr, f_bslt, f_bsmt];   % porosity of structures (nan = do not set)
Tstruct   =       [T_wat, T_sed,  T_plb,      T_mlb, T_imr, T_bslt, T_bsmt];   % temperature of structures (nan = do not set)
Cstruct   =       [C_wat, C_sed,  C_plb,      C_mlb, C_imr, C_bslt, C_bsmt];   % salinity of structures (nan = do not set)

tol     = 1e-8;      % residual tolerance for iterative solver
alpha   = 1.25;      % step size for iterative solver
beta    = 0.99;      % damping parameter for iterative solver

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************
