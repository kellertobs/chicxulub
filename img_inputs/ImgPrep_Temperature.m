%% This is a workbook to create inputs for run_impact.m so that you can have complex temperature distributions
% % Crop image to dimensions you want (differentiate between full crater and half crater)
% % Segment a temperature distribution cross section into units (nUnits_T) using k means
% % Save nUnits binary images where units are True and all other space is False



%% Initial space set up
clc;	% Clear command window.
% clear;	% Delete all variables.
close all;	% Close all figure windows except those created by imtool.
imtool close all;	% Close all figure windows created by imtool.
workspace;	% Make sure the workspace panel is showing.


%% User inputs/options
% run('../usr/par_Sudbury_R_500x500.m');  % Use this parameter file


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
width_T     = width(img);     % get width of original image in pixels
height_T    = height(img);    % get height of original image in pixels

imgCrp    = imcrop(img, [width_T*x_crp, height_T*y_crp, width_T*w_crp, height_T*h_crp]);  % Crop function to select only part of the image defined by rectangle with [xmin ymin width height] REMEMBER: Origin is in top left for MatLab reasons

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
    imgTitle = ['Cluster' num2str(i)];  % Make a title for the image that tells you which cluster it is

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
for l = 1:nUnits_T-1;
    [m,n] = size(c{l});
     for i = 1:m;
        for j = 1:n;
            for p = 1:nUnits_T-l;
            
                if c{l}(i, j) == 1 && (c{l+p}(i, j) == 1);
                    msg = ["replaced" num2str(i) num2str(j)];
%                     disp(msg)
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


%% Assign temperatures to each array 
numRows = nUnits_T;
numCols = 1;

% File to store last values
memoryFile = 'last_temperature_inputs.mat';

% Load previous defaults if available (robust to different saved variable names)
if isfile(memoryFile)
    varsInFile = who('-file', memoryFile);            % list variables saved in the mat
    if ismember('defaultData', varsInFile)
        tmp = load(memoryFile, 'defaultData');
        defaultData = tmp.defaultData;
    elseif ismember('userInput', varsInFile)         % backward compatibility
        tmp = load(memoryFile, 'userInput');
        defaultData = tmp.userInput;
    else
        defaultData = repmat({''}, numRows, 1);
    end
    % Ensure correct size if number of units changed
    if numel(defaultData) ~= numRows
        % If there are fewer stored defaults, pad with blanks; if more, trim
        if numel(defaultData) < numRows
            defaultData = [defaultData; repmat({''}, numRows - numel(defaultData), 1)];
        else
            defaultData = defaultData(1:numRows);
        end
    end
else
    defaultData = repmat({''}, numRows, 1);
end

% Prepare prompts
prompts = cell(numRows, 1);
for r = 1:numRows
    prompts{r} = sprintf('Enter Cluster %d temperature in degrees Celsius:', r);
end

% Show the input dialog box with remembered defaults
userInput = inputdlg(prompts, 'Enter Data', [1 100], defaultData);

% If user canceled the dialog, exit gracefully
if isempty(userInput)
    disp('User canceled input.');
    return;
end

% Save current inputs for next run (save under 'defaultData' for consistency)
defaultData = userInput;   % store as defaultData for next time
save(memoryFile, 'defaultData');

% Parse input into unitData (numRows x numCols cell array)
unitData = cell(numRows, numCols);
for r = 1:numRows
    parts = strsplit(strtrim(userInput{r}), {' ', ','});
    if numel(parts) ~= numCols
        error('Row %d does not contain %d elements.', r, numCols);
    end
    unitData(r, :) = parts;
end


% Create variables in base, display results, compute T_array
T_array = cell(1, numRows);

for r = 1:numRows
    % Convert to numeric
    T_val = str2double(unitData{r});
    if isnan(T_val)
        error('Invalid numeric input at row %d.', r);
    end
    
    % Assign actual numeric value into base workspace variable T_val_r
    assignin('base', sprintf('T_val_%d', r), T_val);
    
    % Print to command window (nice formatting)
    fprintf('T_val_%d = %g\n', r, T_val);
    
    % Compute T_array element (assumes c{r} exists and is compatible)
    T_array{r} = T_val * c{r};
end

% Sum arrays together
T_array2 = 0;
for i = 1:numRows
    T_array2 = T_array2 + T_array{i};
end



% Display the result and confirm satisfaction with the temperature map
happy = false;   % flag for user satisfaction

while ~happy
    % Compute and display the summed temperature array
    T_array2 = 0;
    for i = 1:numRows
        T_array2 = T_array2 + T_array{i};
    end

    figure();
    imagesc(T_array2);
    colorbar;
    title('Summed Temperature Array');
    axis equal tight;
    drawnow;

    % Ask user if they are happy
    choice = questdlg('Are you happy with this temperature distribution?', ...
                      'Confirm Temperature Map', ...
                      'Yes, save and continue', 'No, re-enter temperatures', 'Cancel', ...
                      'Yes, save and continue');

    switch choice
        case 'Yes, save and continue'
            happy = true;  % exit loop and save
        case 'No, re-enter temperatures'
            % Re-run the input dialog for new temperatures
            userInput = inputdlg(prompts, 'Re-enter Temperatures', [1 100], defaultData);

            % If user cancels during re-entry, exit gracefully
            if isempty(userInput)
                disp('User canceled input.');
                return;
            end

            % Save new defaults
            defaultData = userInput;
            save(memoryFile, 'defaultData');

            % Recalculate T_array based on new inputs
            for r = 1:numRows
                T_val = str2double(userInput{r});
                if isnan(T_val)
                    error('Invalid numeric input at row %d.', r);
                end
                assignin('base', sprintf('T_val_%d', r), T_val);
                fprintf('T_val_%d = %g\n', r, T_val);
                T_array{r} = T_val * c{r};
            end

        case 'Cancel'
            disp('User canceled process.');
            return;
    end
end

% === Save temperature array once user is satisfied ===
filename = [foldername '/' projectName '_' num2str(Nx) 'x' num2str(Nz) '_TArray.mat'];    % Specify filename
save(filename, 'T_array2');
fprintf('Temperature array saved to %s\n', filename);

