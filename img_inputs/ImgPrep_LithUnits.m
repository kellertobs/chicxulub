% %% TO DO LIST


%% Initial space set up
clc;	% Clear command window.
clear;	% Delete all variables.
close all;	% Close all figure windows except those created by imtool.
imtool close all;	% Close all figure windows created by imtool.
workspace;	% Make sure the workspace panel is showing.


%% This is a workbook to create inputs for run_impact.m so that you don't need to specify rock units by hand
% % Crop image to dimensions you want (differentiate between full crater and half crater)
% % Segment a geological cross section into units (nUnits) using k means
% % Save nUnits binary images where units are True and all other space is False
% % Save nUnits binary images where units are False and all other space is True ***(WHY? Because you need it sometimes to fill in images that had transparency, for example, water or air above crater)
% % Inputs that you need:
% %     - image of the cross section with different units in different colors and/or textures. Probably best is solid single color fill
% %     - dimensions of the grid size in run_impact.m (N)
% %     - number of expected geological units (need to set for k means segmentation)


%% Initial user inputs/options
imgName = "Lonar/Lonar_Full_Lithology.png"  % Filename of the image of the cross section for input
N = 200                     % dimensions for number of pixels  (should be same as grid size in Run Impact - EXCEPT NEED TO ADD 2 TO MAKE IT MATCH INDSTRUCT BUILDER, NOT SURE WHY)
M = N                       % dimensions for depth of pixels in case this ever stops being square
nUnits = 7;                 % Need to specify number of geological units here; think they're subsequently clustered with 1 being highest area and 5 being lowest
fontSize = 16;              % Fontsize for plotting images while processing

outdir = '../img_inputs/Lonar/'
projectName = 'Lonar_Full_01' % Specify project name so files will saved with some info
foldername = [outdir projectName '_' num2str(N) 'x' num2str(M)];    % Specify foldername for output
mkdir (sprintf(foldername));    % Make the specified directory

smallestAcceptableArea = 50;   % Used for removing tiny bits of white pixels from binary images. Keep areas only if they're bigger than this number in pixels MAY NEED TO ADJUST THIS DEPENDING ON WHAT TYPE OF IMAGE
se = strel('disk', 10);         % Used for removing tiny bits of black pixels from binary images. Make a structuring element to fill in holes in white areas ('shape', size in pixels) MAY NEED TO ADJUST THIS DEPENDING ON WHAT TYPE OF IMAGE


%% Make two figures
% f1 = figure;    % Make a figure for plotting all of the final units
f2 = figure;    % Make a figure for plotting each step of the process
set(gcf, 'Position', get(0, 'ScreenSize')); % Maximize the figure so you can see the steps


%% Import  image
img = imread(imgName);

subplot(3, 4, 1);   % Where to plot the original image in the figure
imshow(img);        % Show the image
title("Original Image RGB");    % Title for the image
drawnow;            % Make it display immediately.


%% Crop the image to whatever size you want
width = width(img);     % get width of original image in pixels
height = height(img);   % get height of original image in pixels

imgCrp = imcrop(img, [0, 0, width, width]);  % Crop function to select only part of the image defined by rectangle with [xmin ymin width height] REMEMBER: Origin is in top left for MatLab reasons

subplot(3, 4, 2);   % Where to plot the cropped image
imshow(imgCrp);     % Show the image
title("Cropped Image RGB");    % Title for the image
drawnow;            % Make it display immediately.


%% Segment image into nUnits (number of geological units specified above)
unitLabels = imsegkmeans(imgCrp, nUnits);          % Segment image
unitOverlay = labeloverlay(imgCrp, unitLabels);    % Overlay segmented images on greyscale image

subplot(3, 4, 3);   % Where to plot the segmented image
imshow(unitOverlay);    % Show the image
title("Segmented Image");   % Title for the image
drawnow;            % Make it display immediately.


%% Split, clean, and save segments into individual binary images
for i = 1:nUnits          
    figure(f2)  % Select Figure 2 for plotting each step in image processing

    % % Split image into different components
    mask = unitLabels == i;             % make a mask for pixels where unitLabels created above match i unit
    imgClstr = imgCrp.*uint8(mask);    % CAN'T REMEMBER WHAT THIS DOES - groups together all pixels from one cluster/segment maybe?
    imgTitle = ['Cluster_' num2str(i)]  % Make a title for the image that tells you which cluster it is

    subplot(3, 4, 4);   % Where to plot the segmented image 
    imshow(imgClstr);   % Show the image
    title(imgTitle);    % Title for the image
    drawnow;            % Make it display immediately.
        

    %% Make greyscale - make greyscale in because you need to 2D array, not 3D array
    imgGrey = im2gray(imgClstr);  % make image greyscale 
    
    subplot(3, 4, 5);   % Where to plot the image
    imshow(imgGrey);    % Show the image
    title("Cropped Image Greyscale");      % Title for the image
    drawnow;            % Make it display immediately


    %% Make cluster a binary image
    imgBW = imbinarize(imgGrey);
    
    subplot(3, 4, 6);   % Where to plot the segmented image      
    imshow(imgBW);      % Show the image
    title("Binary image of cluster");    % Title for the image
    drawnow;            % Make it display immediately.
    

    %% Filter out small small white bits in black
    imgBW2 = (bwareaopen(imgBW, smallestAcceptableArea));

    subplot(3, 4, 7);   % Where to plot the segmented image      
    imshow(imgBW2);     % Show the image
    title("Remove small white bits");    % Title for the image
    drawnow;            % Make it display immediately.
    

    %% Get rid of holes in the white area (remove black bits in white)
    imgBW3 = imclose(imgBW2,se);    % fill in black bits using above made structuring element

    subplot(3, 4, 8);   % Where to plot the segmented image        
    imshow(imgBW3);     % Show the image
    title("Remove small black bits");    % Title for the image
    drawnow;            % Make it display immediately.


    %% Thicken to help with units not matching up once in run_impact
    imgBW4 = bwmorph(imgBW3, 'thicken', 2);     % Add a few pixels to the edges of each unit to fill in gaps between units

    subplot(3,4, 9);    % Where to plot image
    imshow(imgBW4)      % Show the image
    title("thicken");   % Title for the image
    drawnow;            % Make it display immediately


    %% Pixelate
    imgBW5 = imresize(imgBW4, [N M]);   % Pixelate the image to appropriate grid size
    subplot(3, 4, 10);   % Where to plot the segmented image      
    imshow(imgBW5);      % Show the image
    title("Pixelated");  % Title for the image
    drawnow;             % Make it display immediately
    
    d{i} = imgClstr;
    c{i} = imgBW5;  
end

%% Replacing overlapping pixels
for l = 1:nUnits-1
    [m,n] = size(c{l});
     for i = 1:m
        for j = 1:n
            for p = 1:nUnits-l
            
                if c{l}(i, j) == 1 && (c{l+p}(i, j) == 1)
                    msg = ["replaced" num2str(i) num2str(j)]
                    disp(msg)
                    c{l+p}(i,j) = 0;
                end
            end
        end
     end 
end
                

%% Save images
for l = 1:nUnits
    imgTitle = ['Lith_' num2str(l)];  % Make a title for the image that tells you which cluster it is
    filename = [foldername '/' projectName '_' num2str(N) 'x' num2str(M) '_' imgTitle '.png'];    % Specify filename
    imwrite(c{l}, filename);
    msg1 = ['Saved_' imgTitle '.png'];
    disp(msg1)

% %% Invert images and save 
%     imgInv = imcomplement(c{l});    % Invert image in case you need them (depends on what kind of starting image you have)
%     filename = [outdir projectName '_' num2str(N) 'x' num2str(M) '_' imgTitle '_inv.png'];    % Specify filename
%     imwrite(imgInv, filename);         % Save inverted image
%     msg2 = ['Saved_' imgTitle '_inv.png'];
%     disp(msg2)
end


%% Plot each unit/cluster on Figure 1
f1 = figure;    % Make a figure for plotting all of the final units

figure(f1);         % Select Figure 1 to plot each cluster on one figure
subplot(3,4,1);     % Plot cropped original image to the first position on Figure 1
imshow(img)
title("Original RGB");  % Title image

for i = 1:nUnits
    imgTitle = ['Lith ' num2str(i)];  % Make a title for the image that tells you which cluster it is
    subplot(3,4,i+1);     % Each cluster goes in it's own subplot
    imshow(d{i});     % Show image
    title(imgTitle);    % Title for the image
    drawnow
end

filename = [foldername '/' projectName '_' num2str(N) 'x' num2str(M) '_LithClusters.png'];    % Specify filename
saveas(f1, filename)   ;      % Save figure of all clusters
% msg3 = ['Saved_AllClusters.png']
% disp(msg3)
