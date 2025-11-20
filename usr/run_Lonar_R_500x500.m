clear; close all; clc;

par_default;

% SET MODEL PARAMETERS
runID   = 'Sudbury_R'; % run identifier tag
nout    = 20;       % print output every 'nop' steps
svout   = 1;        % save figures and data to file (1)

% set domain parameters
Nz      = 500;      % num. grid size in z-direction
Nx      = 500;      % num. grid size in z-direction 

Nzi     = 100;      % num. grid size in z-direction
Nxi     = 100;      % num. grid size in x-direction
D       = 20e3;      % physical domain depth [m]

indir   = '../img_inputs/Sudbury/Sudbury_R_500x500/'; % input directory for arrays
outdir  = '../out'; % output directory 

%% Make background temperature and porosity from an array (make array using ImgPrep_Temperature.m and ImgPrep_Porosity.m)
TArr   = load([indir 'Sudbury_R_500x500_TArray.mat']);
TArray = TArr.T_array2;

addpath ../src
for i = 1:300
    TArray = TArray + diffus(TArray,ones(size(TArray))/8,1,[1,2],BC_VP);
end

% fArr = load([indir 'Lonar_Right_01_200x200_fArray.mat']);
% fArray = fArr.f_array2;

%% Make lithologic units from binary images (2D arrays) prepared using ImgPrep_LithUnits.m
run('../src/lith_input')

%% set initial condition parameters
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


smth    = 5; % smoothness of initial fields

% set boundary conditions
BC_T    = {{'flux',[0,0.03]},'closed'};
BC_C    = {'closed','closed'};
BC_V    = {'closed','closed'};
BC_VP   = {'closed','closed'};

% downsample image inputs
ho = D./Nz;                       % grid spacing [m]
W  = D*Nx/Nz;                     % domain width [m]
xo = linspace(ho/2,W-ho/2,Nx);    % x-coordinate vector
zo = linspace(ho/2,D-ho/2,Nz);    % z-coordinate vector
[Xo,Zo] = meshgrid(xo,zo);           % coordinate arrays

hi = D./Nzi;                      % grid spacing [m]
W  = D*Nxi/Nzi;                   % domain width [m]
xi = linspace(hi/2,W-hi/2,Nxi);   % x-coordinate vector
zi = linspace(hi/2,D-hi/2,Nzi);   % z-coordinate vector
[Xi,Zi] = meshgrid(xi,zi);           % coordinate arrays

TArray = interp2(Xo,Zo,double(TArray),Xi,Zi);


wat_evolve = false;              % evolve water as well-mixed reservoir; else keep T,C constant
tau_eqlb   = 1*3600*24*365.25;   % water-air thermal equilibration time 

wat   = zeros(Nzi,Nxi);                 % input image has wrong values along topmost row, please fix
air   = zeros(Nzi,Nxi);                 % input image has wrong values along topmost row, please fix

%% --- INITIALIZE STORAGE ---
numLith = height(summaryTable);

unit = [];            % to hold all lith arrays stacked in 3rd dimension
fstruct = zeros(1, numLith);
Tstruct = zeros(1, numLith);
Cstruct = zeros(1, numLith);

%% --- LOOP THROUGH EACH LITHOLOGY IMAGE ---
for j = 1:numLith

    % --- Read lithology image file ---
    lithFile = summaryTable.FileName{j};
    lithPath = fullfile(folderPath, lithFile);

    % Try to find a valid image extension if needed
    if ~isfile(lithPath)
        exts = {'.png','.jpg','.jpeg','.tif','.tiff','.bmp'};
        found = false;
        for e = 1:numel(exts)
            testpath = fullfile(indir, [lithFile exts{e}]);
            if isfile(testpath)
                lithPath = testpath;
                found = true;
                break;
            end
        end
        if ~found
            warning('File "%s" not found. Skipping.', lithFile);
            continue;
        end
    end

    lith_j = imread(lithPath);

    % Convert to double for interpolation
    lith_j = double(lith_j);

    % Downsample / resample to target grid size
    
    targetRows = numel(Zi(:,1));
    targetCols = numel(Xi(1,:));

    lith_j = imresize(lith_j, [targetRows, targetCols], 'bilinear');
    lith_j = round(lith_j);  % optional rounding if binary/integer mask

    % Add as new layer in 3D array
    if isempty(unit)
        unit = lith_j; % initialize
    else
        unit = cat(3, unit, lith_j);
    end

    % --- Assign property values from summary table ---
    fstruct(j) = summaryTable.Porosity_f(j);
    Tstruct(j) = summaryTable.Temperature_T(j);
    Cstruct(j) = summaryTable.Salinity_C(j);

    if strcmp(summaryTable.UnitName{j}, 'Water')
        wat = lith_j;
    end

    if strcmp(summaryTable.UnitName{j}, 'Air')
        air = lith_j;
    end


end

Nx = Nxi;
Nz = Nzi;

%% --- DISPLAY RESULTS ---
fprintf('\n Loaded %d lithology layers.\n', size(unit,3));
fprintf('   Created arrays:\n');
fprintf('   unit   : %d×%d×%d\n', size(unit));
fprintf('   fstruct: [%s]\n', num2str(fstruct));
fprintf('   Tstruct: [%s]\n', num2str(Tstruct));
fprintf('   Cstruct: [%s]\n', num2str(Cstruct));

%% Parameters for iterative solver
tol     = 1e-7;      % residual tolerance for iterative solver
alpha   = 1.25;      % step size for iterative solver
beta    = 0.99;      % damping parameter for iterative solver
maxit   = 5e3;       % maximum number of iterations

%*****  RUN CHICXULUB MODEL  *************************************************
run('../src/main')
%**************************************************************************
