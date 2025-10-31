%% IMAGE PREP SCRIPT
%  For turning image files into usable inputs for hydrothermal model.
%  Requires lithology image and temperature distribution image. 
%  Default is to run everything at 500x500 grid size and downsample in
%  run_Impact.m

%% Initial space set up
clc;	% Clear command window.
clear;	% Delete all variables.
close all;	% Close all figure windows except those created by imtool.
imtool close all;	% Close all figure windows created by imtool.
workspace;	% Make sure the workspace panel is showing.

%% SET MODEL PARAMETERS

runID   = 'Lonar_R'; % run identifier tag

% set domain parameters
Nz      = 400;       % num. grid size in z-direction
Nx      = 400;       % num. grid size in z-direction


%% Set conditions for Img_Prep files
outdir_ImgPrep = '../img_inputs/Lonar/'; 

imgName_Lith   = "../img_inputs/Lonar/Lonar_Full_Lithology_watOnly.png";     % Filename of the lithology distribution
nUnits_Lith    = 7;                                                % Need to specify number of Lithological units here 

imgName_T      = "../img_inputs/Lonar/Lonar_Full_Temperature.png";        % Filename of the temperature distribution  
nUnits_T       = 8;                                                 % Need to specify number of Temperature units here  

%% Crop function to select only part of the image defined by rectangle with [xmin ymin width height] REMEMBER: Origin is in top left for MatLab reasons 
% These are scaling factors for splitting the image as fractions of width and height of the image. 
% Some basic options:
%     - Full image:         (0, 0, 1, 1)
%     - Top right of image: (0.5, 0, 0.5, 0.5)
%     - Top left of image:  (0, 0, 0.5, 0.5)

x_crp = 0.5;  
y_crp = 0;
w_crp = 0.5;
h_crp = 0.5;

%% RUN IMAGE PREP
run('../img_inputs/ImgPrep_LithUnits')
run('../img_inputs/ImgPrep_Temperature')