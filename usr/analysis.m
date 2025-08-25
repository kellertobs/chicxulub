% indir   = '../img_inputs/Lonar/Lonar_Right_05/'; % input directory for arrays

% filename = 10 '.mat'
BaseName1 = '../out/Lonar_Right_04_wide_4/Lonar_Right_04_wide_4_';
BaseName2 = '../out/Lonar_Right_04_wide_5/Lonar_Right_04_wide_5_';

% set domain parameters
Nz      = 100;      % num. grid size in z-direction
Nx      = 200;       % num. grid size in z-direction
D       = 1e3;      % phys. domain depth [m]

plt_depth = 250; % depth from which you want to plot the data (T, etc) - need to be mindful of what will translate to a whole number and/or make a rounding thing so matlab makes a whole number
       
cell_w = D/Nx;
depth  = cell_w.*Nz;
y1 = plt_depth/cell_w;




% Define the folder where the .mat files are located
folder = '../out/Lonar_Right_04_wide_5/'; % Replace with your folder path
filePattern = fullfile(folder, '*.mat'); % Pattern to match .mat files
matFiles = dir(filePattern); % List of .mat files

% Initialize arrays to store time and max values
time = [];
maxValues = [];

% Loop through each .mat file
lngth = length(matFiles)-2
for k=1:lngth
    fullFileName=[BaseName,num2str(k),'.mat']

    fprintf(1, 'Now reading %s\n', fullFileName);
    data = load(fullFileName);
    
    % Assuming the data contains a variable named 'yourVariable' which is a matrix
    % Change 'yourVariable' and 'rowNumber' to match your data structure
    yourVariable = data.T; 
    rowNumber = y1; % Replace with your specific row number
    
    % Extract the maximum value from the specified row
    maxValue = max(yourVariable(rowNumber, :));
    
    % Assume each file name contains the time information as 'timeX.mat'
    % Extract time from the file name (e.g., 'time1.mat' -> 1)
    timeValue = data.time;
    
    % Append the time and max value to the arrays
    time = [time; timeValue];
    maxValues = [maxValues; maxValue];
end

% Plot max value vs time
figure;
plot(time, maxValues, '-', 'LineWidth', 2);

plotDepth = num2str(plt_depth);
plotTitle = sprintf('Max Temperature vs Time at %s m depth', plotDepth);
title(plotTitle);

grid on;


% Save the plot
saveas(gcf, 'MaxValue_vs_Time.png');




% figure;


% for k=1:40
%     FileName=[BaseName,num2str(k),'.mat']
%     
%     Data = load(FileName)
% 
%     T = Data.T  % Which data do you want to plot?
% %     T_max = max(T);
%     T_max = max(T,[],2)
% %     M = max(A,[],dim)
% 
%     time = Data.time;
% 
%     x = time;
% 
%     y = T_max(y1);
% %     y = T_max;
%     
%     plot(x,y);
%     hold on
%     title("T vs time");
%     xlabel("Time [s]");
%     ylabel("T_max [C]");
% %     xlim([0,2000]);
% %     ylim([0,200]);
% end







% for k=1:25
%     FileName=[BaseName,num2str(k),'.mat']
%     
%     Data = load(FileName)
% 
%     T = Data.T  % Which data do you want to plot?
% 
%     x = linspace(0,D,Nx);
%     y = T(y1,:);
%     
%     plot(x,y);
%     hold on
%     title("Temperature vs distance");
%     xlabel("Distance [m]");
%     ylabel("Temperature [C]");
%     xlim([0,2000]);
%     ylim([0,200]);
% end


% Data = load(FilenName)
% 
% T = Data.T  % Which data do you want to plot?
% 
% % set domain parameters
% Nz      = 100;      % num. grid size in z-direction
% Nx      = 200;       % num. grid size in z-direction
% D       = 2e3;      % phys. domain depth [m]
% 
% plt_depth = 500; % depth from which you want to plot the data (T, etc) - need to be mindful of what will translate to a whole number and/or make a rounding thing so matlab makes a whole number
% 
% 
% 
% cell_w = D/Nx;
% depth  = cell_w.*Nz;
% 
% y1 = plt_depth/cell_w;
% 
% figure;
% x = linspace(0,D,Nx);
% y = T(y1,:);
% 
% plot(x,y);
% title("Temperature vs distance");
% xlabel("Distance [m]");
% ylabel("Temperature [C]");
% xlim([0,2000]);
% ylim([0,200]);

