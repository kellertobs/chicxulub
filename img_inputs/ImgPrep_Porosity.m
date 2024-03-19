%% This is a workbook to create inputs for run_impact.m so that you can have complex porosity distributions


%% Initial space set up
clc;	% Clear command window.
clear;	% Delete all variables.
close all;	% Close all figure windows except those created by imtool.
imtool close all;	% Close all figure windows created by imtool.
workspace;	% Make sure the workspace panel is showing.


%% User inputs/options
run('../usr/par_Kardla__Right_01');  % Use this parameter file


% % Adjust these parameters if you're not happy with end pixelated images
% (watching each step plot on the first figure can tell you where you have problems (might be over/undersmoothing things)
smallestAcceptableArea = 50;          % Used for removing tiny bits of white pixels from binary images. Keep areas only if they're bigger than this number in pixels MAY NEED TO ADJUST THIS DEPENDING ON WHAT TYPE OF IMAGE
se        = strel('disk', 10);        % Used for removing tiny bits of black pixels from binary images. Make a structuring element to fill in holes in white areas ('shape', size in pixels) MAY NEED TO ADJUST THIS DEPENDING ON WHAT TYPE OF IMAGE
thck      = 3                         % Used for making units slightly thicker, so that there's not holes          

projectName    = runID         % Specify project name so files will saved with some info
foldername     = [outdir_ImgPrep projectName '_' num2str(N) 'x' num2str(M)];    % Specify foldername for output
mkdir (sprintf(foldername));   % Make the specified directory


%% Make a figure to plot each step of the processing process
f2        = figure;                           % Make a figure for plotting each step of the process
set(gcf, 'Position', get(0, 'ScreenSize'));   % Maximize the figure so you can see the steps


%% Import  image
img       = imread(imgName_f);

subplot(3, 4, 1);   % Where to plot the original image in the figure
imshow(img);        % Show the image
title("Original Image");    % Title for the image
drawnow;            % Make it display immediately.


%% Crop the image to whatever size you want
width     = width(img);     % get width of original image in pixels
height    = height(img);    % get height of original image in pixels

imgCrp    = imcrop(img, [width*x_crp, height*y_crp, width*w_crp, height*h_crp]);  % Crop function to select only part of the image defined by rectangle with [xmin ymin width height] REMEMBER: Origin is in top left for MatLab reasons

subplot(3, 4, 2);   % Where to plot the cropped image
imshow(imgCrp);     % Show the image
title("Cropped Image");    % Title for the image
drawnow;            % Make it display immediately.



%% Segment image into nUnits (number of geological units specified above)
unitLabels = imsegkmeans(imgCrp, nUnits_f);          % Segment image
unitOverlay = labeloverlay(imgCrp, unitLabels);    % Overlay segmented images on greyscale image

subplot(3, 4, 3);   % Where to plot the segmented image
imshow(unitOverlay);    % Show the image
title("Segmented Image");   % Title for the image
drawnow;            % Make it display immediately.


%% Split, clean, and save segments into individual binary images
for i = 1:nUnits_f          
    figure(f2)  % Select Figure 2 for plotting each step in image processing

    % % Split image into different components
    mask = unitLabels == i;             % make a mask for pixels where unitLabels created above match i unit
    imgClstr = imgCrp.*uint8(mask);    % CAN'T REMEMBER WHAT THIS DOES - groups together all pixels from one cluster/segment maybe?
    imgTitle = ['Cluster' num2str(i)]  % Make a title for the image that tells you which cluster it is

    subplot(3, 4, 4);   % Where to plot the segmented image 
    imshow(imgClstr);   % Show the image
    title(imgTitle);    % Title for the image
    drawnow;            % Make it display immediately

        
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
    drawnow;            % Make it display immediately
    

    %% Filter out small small white bits in black
    imgBW2 = (bwareaopen(imgBW, smallestAcceptableArea));

    subplot(3, 4, 7);   % Where to plot the segmented image      
    imshow(imgBW2);     % Show the image
    title("Remove small white bits");    % Title for the image
    drawnow;            % Make it display immediately
    

    %% Get rid of holes in the white area (remove black bits in white)
    imgBW3 = imclose(imgBW2,se);    % fill in black bits using above made structuring element

    subplot(3, 4, 8);   % Where to plot the segmented image        
    imshow(imgBW3);     % Show the image
    title("Remove small black bits");    % Title for the image
    drawnow;            % Make it display immediately


    %% Thicken to help with units not matching up once in run_impact
    imgBW4 = bwmorph(imgBW3, 'thicken', thck);     % Add a few pixels to the edges of each unit to fill in gaps between units

    subplot(3,4, 9);    % Where to plot image
    imshow(imgBW4)      % Show the image
    title("Thicken");   % Title for the image
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
for l = 1:nUnits_f-1
    [m,n] = size(c{l});
     for i = 1:m;
        for j = 1:n;
            for p = 1:nUnits_f-l;
            
                if c{l}(i, j) == 1 && (c{l+p}(i, j) == 1);
                    msg = ["replaced" num2str(i) num2str(j)];
                    disp(msg);
                    c{l+p}(i,j) = 0;
                end
            end
        end
     end 
end

%% Plot each unit/cluster on a new figure
f1 = figure;    % Make a figure for plotting all of the final units
set(gcf, 'Position', get(0, 'ScreenSize'));   % Maximize the figure so you can see the steps

col = ceil((nUnits_f+1)/3);

figure(f1);         % Select Figure 1 to plot each cluster on one figure
subplot(3,col,1);     % Plot cropped original image to the first position on Figure 1
imshow(img)
title("Original RGB");  % Title image
drawnow

for i = 1:nUnits_f
    imgTitle = ['Cluster' num2str(i)];  % Make a title for the image that tells you which cluster it is
    subplot(3,col,i+1);     % Each cluster goes in it's own subplot
    imshow(d{i});     % Show image
    title(imgTitle);    % Title for the image
    drawnow
end

filename = [foldername '/' projectName '_' num2str(N) 'x' num2str(M) '_fClusters.png'];    % Specify filename
saveas(f1, filename)   ;      % Save figure of all clusters


%% Break so you can assign values to clusters
msg4 = ['Assign temperatures to clusters in ImgPrep_Temperature.m; then type dbcont in command window '];
disp(msg4)
keyboard    

%% Assign porosity to each array (need to do this by person brain at the moment, may be a smarter way)
f_array = []

f_array{1} = 0.2*c{8};
f_array{2} = 0.12*c{11};
f_array{3} = 0.05*c{4};
f_array{4} = 0.04*c{3};
f_array{5} = 0.04*c{10};

f_array{6}  = 0.03*c{7};
f_array{7}  = 0.03*c{5};
f_array{8}  = 0.02*c{1};
f_array{9}  = 0.02*c{9};
f_array{10} = 0.01*c{2};

f_array{11}  = 0.5*c{6};   %% Water or air layer

f_array2 = f_array{1} + f_array{2} + f_array{3} + f_array{4} + f_array{5} + f_array{6} + f_array{7} + f_array{8} + f_array{9} + f_array{10} + f_array{11};


%% Save array of porosities

filename = [foldername '/' projectName '_' num2str(N) 'x' num2str(M) '_fArray.mat'];    % Specify filename
save([filename],'f_array2');

