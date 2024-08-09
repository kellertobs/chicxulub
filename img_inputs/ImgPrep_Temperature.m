%% This is a workbook to create inputs for run_impact.m so that you can have complex temperature distributions


%% Initial space set up
clc;	% Clear command window.
clear;	% Delete all variables.
close all;	% Close all figure windows except those created by imtool.
imtool close all;	% Close all figure windows created by imtool.
workspace;	% Make sure the workspace panel is showing.


%% User inputs/options
run('../usr/par_Lonar_Right_04_4.m');  % Use this parameter file


% % Adjust these parameters if you're not happy with end pixelated images
% (watching each step plot on the first figure can tell you where you have problems (might be over/undersmoothing things)
smallestAcceptableArea = 50;          % Used for removing tiny bits of white pixels from binary images. Keep areas only if they're bigger than this number in pixels MAY NEED TO ADJUST THIS DEPENDING ON WHAT TYPE OF IMAGE
se        = strel('disk', 10);        % Used for removing tiny bits of black pixels from binary images. Make a structuring element to fill in holes in white areas ('shape', size in pixels) MAY NEED TO ADJUST THIS DEPENDING ON WHAT TYPE OF IMAGE
thck      = 5;                        % Used for making units slightly thicker, so that there's not holes          

projectName    = runID         % Specify project name so files will saved with some info
foldername     = [outdir_ImgPrep projectName '_' num2str(Nx) 'x' num2str(Nz)];    % Specify foldername for output
mkdir (sprintf(foldername));   % Make the specified directory


%% Make a figure to plot each step of the processing process
f2        = figure;                           % Make a figure for plotting each step of the process
set(gcf, 'Position', get(0, 'ScreenSize'));   % Maximize the figure so you can see the steps


%% Import  image
img       = imread(imgName_T);

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


%% Segment image into nUnits_T (number of geological units specified above)
unitLabels = imsegkmeans(imgCrp, nUnits_T);          % Segment image
unitOverlay = labeloverlay(imgCrp, unitLabels);    % Overlay segmented images on greyscale image

subplot(3, 4, 3);   % Where to plot the segmented image
imshow(unitOverlay);    % Show the image
title("Segmented Image");   % Title for the image
drawnow;            % Make it display immediately.


%% Split, clean, and save segments into individual binary images
for i = 1:nUnits_T          
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
    imgBW5 = imresize(imgBW4, [Nx Nz]);   % Pixelate the image to appropriate grid size
    subplot(3, 4, 10);   % Where to plot the segmented image      
    imshow(imgBW5);      % Show the image
    title("Pixelated");  % Title for the image
    drawnow;             % Make it display immediately
    
    d{i} = imgClstr;
    c{i} = imgBW5;    
end

%% Replacing overlapping pixels
for l = 1:nUnits_T-1
    [m,n] = size(c{l});
     for i = 1:m;
        for j = 1:n;
            for p = 1:nUnits_T-l;
            
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


col = ceil((nUnits_T+1)/3);

figure(f1);         % Select Figure 1 to plot each cluster on one figure
subplot(3,col,1);     % Plot cropped original image to the first position on Figure 1
imshow(img)
title("Original RGB");  % Title image
drawnow

for i = 1:nUnits_T
    imgTitle = ['Cluster' num2str(i)];  % Make a title for the image that tells you which cluster it is
    subplot(3,col,i+1);     % Each cluster goes in it's own subplot
    imshow(d{i});     % Show image
    title(imgTitle);    % Title for the image
    drawnow
end

filename = [foldername '/' projectName '_' num2str(Nx) 'x' num2str(Nz) '_TClusters.png'];    % Specify filename
saveas(f1, filename)   ;      % Save figure of all clusters


%% Break so you can assign values to clusters
msg4 = ['Assign temperatures to clusters in ImgPrep_Temperature.m; then type dbcont in command window '];
disp(msg4)
keyboard    

%% Assign temperatures to each array (need to do this by person brain at the moment, may be a smarter way)
%% Lonar
% T_array = [];T

T_array{1} = 50*c{1};
T_array{2} = 225*c{3};
T_array{3} = T_wat*c{2};
T_array{4} = 25*c{4};
T_array{5} = 300*c{5};

T_array{6}  = 130*c{6};
T_array{7}  = 80*c{7};
T_array{8}  = 375*c{8};
% T_array{9}  = 700*c{9};


T_array2 = T_array{1} + T_array{2} + T_array{3} + T_array{4} + T_array{5} + T_array{6} + T_array{7} + T_array{8};



%% Boltysh
% T_array = [];T
% 
% T_array{1} = 450*c{1};
% T_array{2} = 250*c{2};
% T_array{3} = 30*c{3};
% T_array{4} = 350*c{4};
% T_array{5} = Twat*c{5};
% 
% T_array{6}  = 120*c{6};
% T_array{7}  = 150*c{7};
% T_array{8}  = 550*c{8};
% T_array{9}  = 900*c{9};
% T_array{10}  = 1200*c{10};
% 
% 
% T_array2 = T_array{1} + T_array{2} + T_array{3} + T_array{4} + T_array{5} + T_array{6} + T_array{7} + T_array{8} + T_array{9} + T_array{10};


% T_array{1} = 50*c{1};
% T_array{2} = 100*c{4};
% T_array{3} = 150*c{10};
% T_array{4} = 200*c{2};
% T_array{5} = 250*c{15};
% 
% T_array{6}  = 300*c{7};
% T_array{7}  = 350*c{11};
% T_array{8}  = 400*c{9};
% T_array{9}  = 450*c{14};
% T_array{10} = 500*c{8};
% 
% T_array{11}  = 550*c{12};
% T_array{12}  = 600*c{5};
% T_array{13}  = 650*c{13};
% T_array{14}  = 700*c{6};
% T_array{15}  = 750*c{17};
% 
% T_array{16}  = 800*c{16};
% T_array{17}  = Twat*c{3};
% 
% T_array2 = T_array{1} + T_array{2} + T_array{3} + T_array{4} + T_array{5} + T_array{6} + T_array{7} + T_array{8} + T_array{9} + T_array{10} + T_array{11} + T_array{12} + T_array{13} + T_array{14} + T_array{15} + T_array{16} + T_array{17};

%% Save array of temperatures

filename = [foldername '/' projectName '_' num2str(Nx) 'x' num2str(Nz) '_TArray.mat'];    % Specify filename
save([filename],'T_array2')

