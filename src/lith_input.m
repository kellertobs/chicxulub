folderPath = sprintf([indir 'LithClusters']); % input directory for arrays
indir = folderPath;
 
%% Check for existing summary table
existingSummaryFile = dir(fullfile(folderPath, [runID '_inputsummary.csv']));

useExisting = false;
previousData = [];

if ~isempty(existingSummaryFile)
    
    summaryTable = readtable(fullfile(folderPath, existingSummaryFile(1).name));
    disp(summaryTable);
    
    % Ask user what to do
    choice = questdlg('An existing summary table was found. What would you like to do?', ...
        'Existing Summary Table Found', ...
        'Use Existing Data','Enter/Update Data','Cancel','Enter/Update Data');

    switch choice
        case 'Use Existing Data'    % Read existing summary
            useExisting = true;

        case 'Enter/Update Data'    % Load existing summary for pre-filling
            previousData = readtable(fullfile(folderPath, existingSummaryFile(1).name));

        case 'Cancel'               % Cancel
            error('Operation cancelled by user.');
    end
end

%% If using existing data, skip input collection
if useExisting
    disp('Using existing summary table:');
    disp(summaryTable);
    return;
end

%% If not using existing data, or no data exists:
% Find all images in folder
imgExts = {'*.jpg','*.jpeg','*.png','*.tif','*.tiff','*.bmp'};  % a list of all files extensions that might be images
fileList = [];                                                  % make an empty array to contain filenames
for k = 1:numel(imgExts)                                        % for k = 1 to k = the number of elements in imgExts list
    fileList = [fileList; dir(fullfile(folderPath, imgExts{k}))];   % add to fileList by getting a list of all files in the selected folder that have image extensions
end

if isempty(fileList)
    error('No image files found in the selected folder.');
end

% Identify all base (_Original) images
baseFiles = fileList(contains({fileList.name}, '_Original'));   
if isempty(baseFiles)
    error('No "_Original" images found in the selected folder.');
end

summaryData = {};                                               % create an empty cell array to hold results of data input

% Loop through base images
for b = 1:numel(baseFiles)
    [~, baseName, baseExt] = fileparts(baseFiles(b).name);      % split the filename into its name and extension, but ignore the folder
    basePath = fullfile(folderPath, [baseName baseExt]);        

    % Derive common prefix (everything before '_Original')
    prefix = extractBefore(baseName, '_Original');
    if isempty(prefix)
        prefix = baseName;
    end

    % Find overlay images that share this prefix but contain '_Lith_'
    overlayCandidates = fileList(contains({fileList.name}, prefix) & contains({fileList.name}, '_Lith_'));

    % Read base image
    baseImg = imread(basePath);

    % Process each overlay image for this base
    for o = 1:numel(overlayCandidates)
        [~, overlayName, overlayExt] = fileparts(overlayCandidates(o).name);
        overlayPath = fullfile(folderPath, [overlayName overlayExt]);

        % Read overlay image
        overlayImg = imread(overlayPath);
        % Resize overlay if necessary
        if ~isequal(size(baseImg,1),size(overlayImg,1)) || ~isequal(size(baseImg,2),size(overlayImg,2))
            overlayImg = imresize(overlayImg, [size(baseImg,1), size(baseImg,2)]);
        end

        % Create figure
        fig = figure('Name',sprintf('Base: %s | Overlay: %s', baseFiles(b).name, overlayCandidates(o).name), ...
                     'NumberTitle','off', ...
                     'Color','w', ...
                     'Position',[200 200 900 700]);
        ax = axes('Parent',fig);
        imshow(baseImg,'Parent',ax,'InitialMagnification','fit');
        hold(ax,'on');
        hOverlay = imshow(~overlayImg,'Parent',ax);
        hOverlay.AlphaData = 0.65;                      %% Adjust image transparency here. 0 = totally transparent; 1 = totally opaque. 0.65 is pretty good
        title(ax,sprintf('%s (Overlayed on %s)', overlayCandidates(o).name, baseFiles(b).name), 'Interpreter','none');

        % Determine default inputs if previous data exists
        defInput = {'', '', '', ''};
        if ~isempty(previousData)
            idx = find(strcmp(previousData.FileName, overlayName));     % Find the index (row) of the existing summary table that matches the current image overlay
            if ~isempty(idx)
                defInput = { ...
                    char(previousData.UnitName(idx)), ...
                    num2str(previousData.Porosity_f(idx)), ...
                    num2str(previousData.Temperature_T(idx)), ...
                    num2str(previousData.Salinity_C(idx)) };
            end
        end

        % Input dialog for this overlay
        prompt = { ...
            sprintf('Enter Unit Name for "%s":', overlayName), ...
            'Enter Porosity (f) [volume fraction, 0.0001 to 1]:', ...
            'Enter Temperature (T) [degrees C]. If using T array enter "NaN":', ...
            'Enter Salinity (C) [weight fraction - sea water = 0.035 (3.5% of mass is salt)]:'};
        dlgTitle = sprintf('Enter/Confirm Data for "%s"', overlayCandidates(o).name);
        dims = [1 60];              

        answer = inputdlg(prompt, dlgTitle, dims, defInput);
        
        if isempty(answer)
            close(fig);
            error('Input cancelled for file "%s".', overlayName);
        end

        % Parse user input
        unitName = strtrim(answer{1});
        f_val = str2double(answer{2});
        T_val = str2double(answer{3});
        C_val = str2double(answer{4});

        if isempty(unitName)
            close(fig);
            error('Unit name missing for file "%s".', overlayName);
        end

        % --- Assign workspace variables ---
        fName = sprintf('f_%s', overlayName);
        TName = sprintf('T_%s', overlayName);
        CName = sprintf('C_%s', overlayName);
        assignin('base', fName, f_val);
        assignin('base', TName, T_val);
        assignin('base', CName, C_val);

        % --- Store summary data ---
        summaryData(end+1,:) = {unitName, overlayName, f_val, T_val, C_val};

        % Close figure
        close(fig);
    end
end

% --- Create summary table ---
summaryTable = cell2table(summaryData, 'VariableNames', {'UnitName','FileName','Porosity_f','Temperature_T','Salinity_C'});

% --- Display and save summary ---
% disp(' ');
title = sprintf('===== %s SUMMARY TABLE =====', runID);
disp(title)
disp(summaryTable);
disp('=========================');


csvPath = fullfile(folderPath, [runID '_inputsummary.csv']);
writetable(summaryTable, csvPath);

fprintf('\nSummary table saved to:\n%s\n', csvPath);
