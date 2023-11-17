% %% TO DO LIST
% - 
% - 


%% Initial space set up
clc;	% Clear command window.
clear;	% Delete all variables.
close all;	% Close all figure windows except those created by imtool.
imtool close all;	% Close all figure windows created by imtool.
workspace;	% Make sure the workspace panel is showing.

%% This is a workbook to create inputs for run_impact.m so that you can have complex temperature distributions
% % Crop image to dimensions you want (differentiate between full crater and half crater)
% % Inputs:
% %     - greyscale image of the temperature distributions, brightest area should be highest temperature, darkest areas should be lowest T

% % Output: N by M array with temperature distributions


%% Initial user inputs/options
imgName = "TestTemp2.png"  % Filename of the image of the cross section

N = 202         % dimensions for number of pixels  (should be same as grid size in Run Impact - EXCEPT NEED TO ADD 2 TO MAKE IT MATCH INDSTRUCT BUILDER, NOT SURE WHY)
M = N           % dimensions for depth of pixels in case this ever stops being square
fontSize = 16;  % Fontsize for plotting images while processing

projectName = 'Kardla_TopRight' % Specify project name so files will saved with some info


%% Make a figure
f1 = figure;        % Make a figure for displaying steps
set(gcf, 'Position', get(0, 'ScreenSize')); % Maximize the figure so you can see the steps


%% Import  image
img = imread(imgName);

subplot(3, 4, 1);   % Where to plot the original image in the figure
imshow(img);        % Show the image
title("Original Image");    % Title for the image
drawnow;            % Make it display immediately.


% %% Crop the image to whatever size you want
% width = width(img);     % get width of original image in pixels
% height = height(img);   % get height of original image in pixels
% 
% img = imcrop(img, [2, 0, 4, 4]);  % Crop function to select only part of the image defined by rectangle with [xmin ymin width height] REMEMBER: Origin is in top left for MatLab reasons
% 
% subplot(3, 4, 2);   % Where to plot the cropped image
% imshow(img);     % Show the image
% title("Cropped Image");    % Title for the image
% drawnow;            % Make it display immediately.


%% Make greyscale - make greyscale in here in case you saved wrong file type when making the input, need to be sure it's 2D array, not 3D array
imgGrey = im2gray(img);  % make image greyscale 

subplot(3, 4, 3);   % Where to plot the image
imshow(imgGrey);    % Show the image
title("Cropped Image Greyscale");      % Title for the image
drawnow;            % Make it display immediately.


% %% Invert image
% imgInv = imcomplement(imgGrey);    % Invert image in case you need them (depends on what kind of starting image you have)
% subplot(3, 4, 4);   % Where to plot the cropped image
% imshow(imgInv);     % Show the image
% title("InvertedGreyScale");      % Title for the image
% drawnow;            % Make it display immediately.


%% Pixelate
imgBW5 = imresize(imgGrey, [N M]);
subplot(3, 4, 4);  % Where to plot the segmented image      
imshow(imgBW5);     % Show the image
title("Pixelated");     % Title for the image
drawnow;            % Make it display immediately.


%% Change all values in array with temperatures

max_grey = max(imgBW5, [], 'all' )  % Maximum luminosity in greyscale image (colour of max T contours)
min_grey = min(imgBW5, [], 'all' )  % Minimum luminosity in greyscale image (colour of min T contours)

max_T = 800    % Maximum Temperature 
min_T = 50      % Minimum Temperature 

n_T   = length(unique(imgBW5))  % Total number of unique values in array

contour_width = (max_T-min_T)/n_T

scaling_factor = double(max_T/max_grey);
TempArray = double(imgBW5);

T_array = TempArray.*scaling_factor;

max_T_array = max(T_array, [], 'all')
min_T_array = min(T_array, [], "all")


%% Save array of temperatures
% save('TestTArray.mat', 'T_array');

runID   = 'impact_2'; % run identifier tag
outdir  = '../img_inputs'; % output directory 

save([outdir,'/', 'TestTAray.mat'],'T_array');

