%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjustable program that produces images (and supplementary plots, if
% relevant) of ultrasonic data.
%
% Parallel computing has been used to speed up the data processing. 
% If one wishes to remove this feature, replace all instances of
% 'parfor' with 'for'.
%
% Inputs
% 1. String of folder paths. The folder should contain ADC .m data files for
%    for each angle of receiver/transmitter data. There should be
%    corresponding calibration folder paths.
% 2. Global setup and image processing parameters.
%
% Outputs
% 1. Ultrasonic image of the scene with transmission data and automated
%    attempts at finding local maxima (anomalies).
% 2. Plot of x-axis data at y-value of highest prominence local maximum.
% 3. Plot of y-axis data at x-value of highest prominence local maximum.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean everything up.
close all; clear;

%%% FOLDER PATH INPUTS %%%
folderPath = './11-14 30air/11-14 air air 1 15sweep 30v/rx/';
folderPathCali = './11-14 30air/11-14 air air 11 15sweep 30v/rx/';
folderPathTx = './11-14 30air/11-14 air air 2 15sweep 30v/';
folderPathTxCali = './11-14 30air/11-14 air air 11 15sweep 30v/';

%%% IMAGE PROCESSING PARAMETERS %%%
function G = dspglobals
    G.INTERPOLATION_TYPE = 'nearest'; % 'nearest', 'natural' for img interp.
    G.VOLTAGE_AUTORANGE = false;       % If false, set the two params below.
    G.MAX_VOLTAGE = 475;              % mV
    G.MIN_VOLTAGE = 20;
    G.TRANSMISSION_VOLTAGE_AUTORANGE = false; % If false, set the two params below.
    G.TRANSMISSION_RECEIVERS_MAX_VOLTAGE = 405;
    G.TRANSMISSION_RECEIVERS_MIN_VOLTAGE = 0;
    G.PROMINENCE_THRESHOLD_MULTIPLIER = 0.95; % Strength of local maxima filters. 0.0 - 1.0.
    G.VOLTAGE_THRESHOLD_MULTIPLIER = 0.2;    % Filter for small reflections. 0.0 - 1.0.
    G.DO_X_MASK = false;      % Whether to filter around each local maxima in the x direction.
    G.YMASK_TOLERANCE = 3;   % How narrowly to filter around each local maxima (y).
    G.XMASK_TOLERANCE = 30;  % How narrowly to filter around each local maxima (x).
end

APPLY_TRANSMISSIONMODE = true;

%%% HARDWARE PARAMETERS %%%
samplingRate = 240385;          % XADC sampling rate in Hz.
pulseFreq = 40000;              % Frequency of the pulse.
angleStepSize = 1;              % Angle step size.
propSpeed = 344;                % Speed of sound in meters/second.
lambda = propSpeed / pulseFreq; % Wavelength.
k = 360 / lambda;
receiverDistance = 0.015;       % Distance between receivers in meters.


%% Setup Code %%

% Calculate interpolation and downsampling factors
interpAngle = angleStepSize/1;
del_t = receiverDistance * sind(interpAngle) / propSpeed;
desired_ratio = (1/samplingRate) / del_t;
[N, D] = rat(desired_ratio, 1e-3);
interpNum = N;
interpDen = D;

% Calculate new sampling rate after interpolation and downsampling
interpolatedSamplingRate = (samplingRate * interpNum) / interpDen;

files = dir(fullfile(folderPath, '*.mat')); % Which files to grab.
filesCali = dir(fullfile(folderPathCali, '*.mat'));
filesTx = struct();
filesTxCali = struct();

% Transmission mode folders need special handling.
if (APPLY_TRANSMISSIONMODE)
   controlSequences = {'ctrl_00', 'ctrl_01', 'ctrl_10', 'ctrl_11'};
    for i = 1:length(controlSequences)
        controlFolder = fullfile(folderPathTx, controlSequences{i});
        controlFolderCali = fullfile(folderPathTxCali, controlSequences{i});
        
        % Check if the control folder exists
        if isfolder(controlFolder) && isfolder(controlFolderCali)
            % Retrieve .mat files
            txFiles = dir(fullfile(controlFolder, '*.mat'));
            txFilesCali = dir(fullfile(controlFolderCali, '*.mat'));
            
            % Store the struct array directly under each control sequence in filesTx
            filesTx.(controlSequences{i}) = txFiles;
            filesTxCali.(controlSequences{i}) = txFilesCali;
        else
            disp(['Folder not found: ', controlFolder, controlFolderCali]);
        end
    end 
    filesTx = sort_files_within_control(filesTx);
    filesTxCali = sort_files_within_control(filesTxCali);
end


process_all_files(APPLY_TRANSMISSIONMODE, true, false, true, folderPath, folderPathCali, folderPathTx, folderPathTxCali, files, filesCali, filesTx, filesTxCali, interpNum, interpDen, samplingRate, interpolatedSamplingRate, propSpeed, receiverDistance);


%% File Processing %%

function process_all_files(doTransmission, doSubtraction, doFilter, doInterpolation, ...
    folderPath, folderPathCali, folderPathTx, folderPathTxCali, files, filesCali, filesTx, filesTxCali, ...
    interpNum, interpDen, samplingRate, interpolatedSamplingRate, speedOfSound, receiverDistance)

    % Extract angles from filenames and sort files by angle
    angles = zeros(size(files));
    for i = 1:length(files)
        filename = files(i).name;
        angleStr = regexp(filename, '-?\d+', 'match');
        angles(i) = str2double(angleStr{end});  % Assumes the angle is the last number in the filename
    end
    [angles, sortedIndices] = sort(angles);
    sortedFiles = files(sortedIndices);
    sortedFilesCali = filesCali(sortedIndices); % Assumes a naming system

    % For plotting
    % Define offset & time step size
    timeStepSize = 1/interpolatedSamplingRate; % Time step in seconds
    
    delayedSigTotal = [];
    delayedSigTotalCali = [];

    % % Bandpass filter design
    lowCutoff = 39000;  % Lower cutoff frequency in Hz
    highCutoff = 41000;  % Upper cutoff frequency in Hz
    order = 4;  % Filter order
    [b, a] = butter(order, [lowCutoff, highCutoff] / (samplingRate / 2), 'bandpass');
    
    % Transmission mode handling
    if(doTransmission)
        receiverAngles = [11, 9, 8, 6, 5, 3, 1, 0, -1, -3, -5, -6, -8, -9, -11, -12];
        maxValues = zeros(16, 1);
        maxValuesCali = zeros(16, 1);
        maxValuesRaw = zeros(16, 1); % deep copy to keep track of the original max values
        controlSequences = fieldnames(filesTx);
        for i = 1:length(controlSequences)
            ctrl = controlSequences{i};
            sortedTxFiles = filesTx.(ctrl);
            sortedTxFilesCali = filesTxCali.(ctrl);
    
            % Map receivers in this control sequence to global receiver indices
            switch ctrl
                case 'ctrl_00'
                    receiverIndices = [1, 2, 3, 4];
                case 'ctrl_01'
                    receiverIndices = [5, 6, 7, 8];
                case 'ctrl_10'
                    receiverIndices = [9, 10, 11, 12];
                case 'ctrl_11'
                    receiverIndices = [13, 14, 15, 16];
            end
    
            numAngles = length(sortedTxFiles);
            numReceiversPerCtrl = 4;
    
            % Initialize arrays to store max data and angles for this control sequence
            maxDataCtrl = zeros(numReceiversPerCtrl, numAngles);
            maxDataCtrlCali = zeros(numReceiversPerCtrl, numAngles);
            anglesCtrl = zeros(1, numAngles);
    
            parfor j = 1:numAngles
                filename = sortedTxFiles(j).name;
                filenameCali = sortedTxFilesCali(j).name;
                filepath = fullfile(sortedTxFiles(j).folder, filename);
                filepathCali = fullfile(sortedTxFilesCali(j).folder, filenameCali);
    
                % Load and process the file
                disp(['Processing ', filepath, ' for control ', ctrl]);
                disp(['    With ', filepathCali]);
    
                % Measurement data
                [r4, r3, r2, r1, ~, angle] = get_receiver_data(filepath, samplingRate);
                r1 = calculate_segmented_max(r1, b, a, samplingRate);
                r2 = calculate_segmented_max(r2, b, a, samplingRate);
                r3 = calculate_segmented_max(r3, b, a, samplingRate);
                r4 = calculate_segmented_max(r4, b, a, samplingRate);
                
                % Calibration data
                [r4c, r3c, r2c, r1c, ~, ~] = get_receiver_data(filepathCali, samplingRate);
                r1c = calculate_segmented_max(r1c, b, a, samplingRate);
                r2c = calculate_segmented_max(r2c, b, a, samplingRate);
                r3c = calculate_segmented_max(r3c, b, a, samplingRate);
                r4c = calculate_segmented_max(r4c, b, a, samplingRate);
    
                % Store the angle
                anglesCtrl(j) = angle;
    
                % Store max values for each receiver in this control sequence
                maxDataCtrl(:, j) = [r1; r2; r3; r4];
                maxDataCtrlCali(:, j) = [r1c; r2c; r3c; r4c];
            end
    
            % Get expected angles for receivers in this control sequence
            expectedAngles = receiverAngles(receiverIndices);
    
            % Initialize array to store max values per control sequence
            maxValuesPerCtrl = zeros(numReceiversPerCtrl, 1);
            maxValuesCaliPerCtrl = zeros(numReceiversPerCtrl, 1);
            maxValuesRawPerCtrl = zeros(numReceiversPerCtrl, 1);
    
            % For each receiver in this control sequence
            for r = 1:numReceiversPerCtrl
                expectedAngle = expectedAngles(r);
    
                % Find indices for expected angle and its neighbors (+/-1 degree)
                idx = find(anglesCtrl == expectedAngle);
                idxMinus1 = find(anglesCtrl == expectedAngle - 1);
                idxMinus2 = find(anglesCtrl == expectedAngle - 2);
                idxPlus1 = find(anglesCtrl == expectedAngle + 1);
                idxPlus2 = find(anglesCtrl == expectedAngle + 2);
                idxs = [idxMinus2, idxMinus1, idx, idxPlus1, idxPlus2];
                idxs = idxs(idxs > 0); % Remove zero or negative indices

                % Collect data and compute the maximum
                if ~isempty(idxs)
                    data = maxDataCtrl(r, idxs);
                    dataCali = maxDataCtrlCali(r, idxs);
                    maxValuesPerCtrl(r) = max((max(dataCali) - max(data)), 0);
                    maxValuesCaliPerCtrl(r) = max(dataCali);
                    maxValuesRawPerCtrl(r) = max(data);
                else
                    maxValuesPerCtrl(r) = NaN; % Handle missing data appropriately
                    maxValuesRawPerCtrl(r) = NaN;
                end
            end
    
            % Store the max values into the global maxValues array
            maxValues(receiverIndices) = maxValuesPerCtrl;
            maxValuesCali(receiverIndices) = maxValuesCaliPerCtrl;
            maxValuesRaw(receiverIndices) = maxValuesRawPerCtrl;
        end
    end

    % This can help us normalize the outer angles
    calibrationDivision = maxValuesCali ./ max(maxValuesCali);
    maxValues = (maxValues ./ calibrationDivision);

    % For shifting the signal back to match the physical location
    timeShift_s = 583e-6;  % 583 microseconds
    %timeShift_s = 0;
    numSamplesShift = round(timeShift_s * interpolatedSamplingRate);
    
    % Process each file in sorted order
    parfor i = 1:length(sortedFiles)
        % File name for this angle
        dataFileName = fullfile(folderPath, sortedFiles(i).name);
        fprintf('Processing %s:\n', dataFileName);

        % Get the data
        [r1, r2, r3, r4, burstSem, angle] = get_receiver_data(dataFileName, samplingRate);

        % Apply the bandpass filter to the data
        r1 = filtfilt(b, a, r1);
        r2 = filtfilt(b, a, r2);
        r3 = filtfilt(b, a, r3);
        r4 = filtfilt(b, a, r4);

        % Average the bursts
        [av1, av2, av3, av4] = average_bursts(r1, r2, r3, r4, burstSem, samplingRate, dataFileName);

        % Calculate necessary shift factor for the given angle
        del_t = receiverDistance * sind(angle) / speedOfSound;
        shiftFactor = abs(del_t / timeStepSize);

        % Apply delays
        [delayedSig, av1, av2, av3, av4] = apply_delays(av1, av2, av3, av4, angle, interpNum, interpDen, shiftFactor);
        
        if(doSubtraction)
            dataFileNameCali = fullfile(folderPathCali, sortedFilesCali(i).name);
            fprintf('     With calibration %s:\n', dataFileNameCali)
            [r1c, r2c, r3c, r4c, burstSemc, ~] = get_receiver_data(dataFileNameCali, samplingRate);
            r1c = filtfilt(b, a, r1c);
            r2c = filtfilt(b, a, r2c);
            r3c = filtfilt(b, a, r3c);
            r4c = filtfilt(b, a, r4c);
            [av1c, av2c, av3c, av4c] = average_bursts(r1c, r2c, r3c, r4c, burstSemc, samplingRate, dataFileNameCali);
            [delayedSigCali, ~, ~, ~, ~] = apply_delays(av1c, av2c, av3c, av4c, angle, interpNum, interpDen, shiftFactor);
            delayedSig = max(delayedSig - delayedSigCali, 0); % Subtraction done here.
        end
    
        % time shift to account for transfer function of transmitters.
        % aligns peak with physical location.
        shiftedSig = zeros(size(delayedSig));
        oldIdx = (1 + numSamplesShift) : length(delayedSig);
        newIdx = 1 : (length(delayedSig) - numSamplesShift);
        shiftedSig(newIdx) = delayedSig(oldIdx);

        delayedSigTotal = [delayedSigTotal, shiftedSig];
    end

    plot_ultrasound_image(delayedSigTotal, angles, timeStepSize, speedOfSound, ...
                          doTransmission, maxValues, maxValuesRaw, receiverAngles);
end


%% Image Generation %%

function plot_ultrasound_image(delayedSigTotal, angles, timeStepSize, speedOfSound, doTransmission, ...
                               maxValues, maxValuesRaw, receiverAngles)
    %% Polar to Cartesian
    % Converts polar coordinates to Cartesian coordinates and plots the ultrasound image.
    [numSamples, numAngles] = size(delayedSigTotal);

    % Define the range for angles and distances.
    maxDistance = (numSamples - 1) * timeStepSize * speedOfSound/2;

    % Initializing the matrices used to store Cartesian data.
    xResolution = 4001;
    yResolution = 2000; % x-axis covers twice the range of the y-axis
    x = linspace(-maxDistance, maxDistance, xResolution) * 100; % Multiplying by 100 to convert to cm.
    y = linspace(0, maxDistance, yResolution) * 100;

    intensity = zeros(yResolution, xResolution);

    edgex1 = [];
    edgex2 = [];
    edgey1 = [];
    edgey2 = [];

    % Map polar data to Cartesian grid and calculate intensity
    for col = 1:numAngles
        theta = angles(col) + 90; % Image is rotated 90 degrees relative to the unit circle.
        for row = 1:numSamples
            distance = row * timeStepSize * speedOfSound/2; 
            xi = distance * cosd(theta);
            yi = distance * sind(theta);

            % Calculate intensity as the absolute difference from the DC offset
            currentIntensity = abs(delayedSigTotal(row, col) - 0);

            % Convert Cartesian coordinates to matrix indices
            xIdx = round((xi + maxDistance) / (2 * maxDistance) * (xResolution - 1)) + 1;
            yIdx = round(yi / maxDistance * (yResolution - 1)) + 1;

            if(col == 1)
                edgex1 = [edgex1, xIdx];
                edgey1 = [edgey1, yIdx];
            elseif(col == numAngles)
                edgex2 = [edgex2, xIdx];
                edgey2 = [edgey2, yIdx];
            end

            % Set intensity
            intensity(yIdx, xIdx) = currentIntensity;
        end
    end

    % Mask out reflection data beyond y = 60 cm
    % As this is beyond where the transmission mode PCB is
    reflectionMask = (y <= 60 - 15); % -15 accounts for the time-shift correction
    intensity(~reflectionMask, :) = 0;

    G = dspglobals();
    
    % Voltage conversion, find maxes/mins
    voltage = double(intensity);
    if (G.VOLTAGE_AUTORANGE)
        maxVoltage = max(voltage(:));
        minVoltage = 0;
    else
        maxVoltage = G.MAX_VOLTAGE;
        minVoltage = G.MIN_VOLTAGE;
    end

    % Interpolation
    [maskY, maskX, mask_intensity] = find(voltage);
    [X, Y] = meshgrid(1:size(voltage, 2), 1:size(voltage, 1));
    intensity_filled = griddata(maskX, maskY, mask_intensity, X, Y, G.INTERPOLATION_TYPE);
    intensity_rayleigh = intensity_filled;
    intensity_filled = uint8(255 * mat2gray(intensity_filled,[minVoltage maxVoltage]));

    for i = 1:length(edgey1)
        intensity_filled(edgey1(i),1:edgex2(i)) = 0;
        intensity_filled(edgey1(i),edgex1(i):size(intensity_filled, 2)) = 0;
    end

    % Apply mask
    intensity_filled(~reflectionMask, :) = 0;
    
    %% Local max detection
    [lmax, prom] = islocalmax2(intensity_filled);
    
    % Filtering lmax based on prominence & intensity in a neighborhood
    promThreshold = G.PROMINENCE_THRESHOLD_MULTIPLIER * max(prom(:));
    peakIntensityThreshold = G.VOLTAGE_THRESHOLD_MULTIPLIER * maxVoltage;
    
    % Define a square neighborhood radius (e.g., 1 means a 3x3 window)
    windowRadius = 10;
    
    % Iterate through each detected local maximum and check its neighborhood
    [localRows, localCols] = find(lmax);
    for i = 1:length(localRows)
        r = localRows(i);
        c = localCols(i);
        
        % Define the boundaries of the neighborhood ensuring they don't go out-of-bounds
        rowStart = max(r - windowRadius, 1);
        rowEnd   = min(r + windowRadius, size(intensity_filled, 1));
        colStart = max(c - windowRadius, 1);
        colEnd   = min(c + windowRadius, size(intensity_filled, 2));
        
        % Extract the neighborhood patch
        neighborhood = voltage(rowStart:rowEnd, colStart:colEnd);
        
        % Reject the maximum if any value in the patch falls below the threshold
        if any(neighborhood(:) >= peakIntensityThreshold)
            lmax(r, c) = true;
        else
            lmax(r, c) = false;
        end
    end

    % Now, also filter based on the prominence threshold.
    lmax = lmax & (prom >= promThreshold);
    [xGrid, yGrid] = meshgrid(1:size(intensity_filled, 2), 1:size(intensity_filled, 1));
    [lmaxY, lmaxX] = find(lmax); % Get row (y) and column (x) indices of lmax points
    toleranceY = G.YMASK_TOLERANCE;
    toleranceX = G.XMASK_TOLERANCE;

    radiusMask = false(size(intensity_filled));

    % Iterate over each local maximum and apply a range around it
    for i = 1:length(lmaxY)
        currentY= lmaxY(i);
        currentX = lmaxX(i);
        % Create a y-based (and optionally x-based) mask for this point
        if (G.DO_X_MASK)
            currentMask = (yGrid >= (currentY - toleranceY)) & (yGrid <= (currentY + toleranceY)) & (xGrid >= (currentX - toleranceX)) & (xGrid <= (currentX + toleranceX));
        else
            currentMask = (yGrid >= (currentY - toleranceY)) & (yGrid <= (currentY + toleranceY));
        end
        radiusMask = radiusMask | currentMask;
    end

    % Apply the mask to the intensity matrix
    filtered_intensity = zeros(size(intensity_filled));
    filtered_intensity(radiusMask) = intensity_filled(radiusMask);
    intensity_filled = filtered_intensity;

%% Plotting
    if (doTransmission)
        figure('Color', 'w');
        
        imagesc(x, y, intensity_filled);
        colormap(gray);
        hcb1 = colorbar('eastoutside');
        caxis([0 255]);
        
        set(hcb1, 'Ticks', [0, 255], 'TickLabels', {'0.0', '1.0'}, ...
                  'FontSize', 14, 'Color', 'k');
        
        ylabel(hcb1, 'Reflection Intensity (Normalized)', 'FontSize', 16, ...
               'FontWeight', 'bold', 'Color', 'k');
        
        set(gca, 'YDir', 'normal');
        
        xlabel('Distance (cm)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
        ylabel('Distance (cm)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
        
        set(gca, 'FontSize', 14, 'FontName', 'Arial', 'LineWidth', 1.5, ...
                 'XColor', 'k', 'YColor', 'k', 'TickDir', 'out');
        box on; 
        
        hold on;
        % Overlay red dots for local maxima
        plot(x(lmaxX), y(lmaxY), 'r.', 'MarkerSize', 10);
        
        % Use the known receiver angles for plotting
        receiverAnglesForPlotting = [10.9, 9.5, 8, 6.6, 5.1, 3.7, 2.2, 0.7, ...
                                    -0.7, -2.2, -3.7, -5.1, -6.6, -8, -9.5, -10.9];
        
        if (G.TRANSMISSION_VOLTAGE_AUTORANGE)
            trans_max = max(maxValues);
            trans_min = 0;
        else
            trans_min = G.TRANSMISSION_RECEIVERS_MIN_VOLTAGE;
            trans_max = G.TRANSMISSION_RECEIVERS_MAX_VOLTAGE;
        end
        
        y_transmission = 65; % in cm, set to 65 for our setup
        numReceivers = length(receiverAnglesForPlotting);
        transmission_x = zeros(numReceivers, 1);
        transmission_y = y_transmission * ones(numReceivers, 1); % All at y_transmission
        transmission_values = maxValues(:); % Ensure it's a column vector
        transmission_values_clipped = min(max(transmission_values, trans_min), trans_max);
        
        % Compute x positions from receiver angles
        for idx = 1:numReceivers
            receiverAngle = receiverAnglesForPlotting(idx); % Use plotting angles
            % Convert angle to x position at y_transmission
            transmission_x(idx) = y_transmission * tand(receiverAngle);
        end
            
        norm_values = (transmission_values_clipped - trans_min) / (trans_max - trans_min);
        
        % Handle potential division by zero if all values are the same
        if any(isnan(norm_values))
            norm_values = zeros(size(transmission_values)); % Assign zero to all if normalization fails
        end
        
        cmap = jet(256);
        norm_values_idx = round(norm_values * (size(cmap, 1) - 1)) + 1;
        
        norm_values_idx = max(1, min(size(cmap, 1), norm_values_idx));
        
        norm_values_idx = size(cmap, 1) - norm_values_idx + 1;

        colors = cmap(norm_values_idx, :);
        
        scatter(transmission_x, transmission_y, 150, colors, 'filled', 'MarkerEdgeColor', 'k');
        
        xlim([-18, 18]); 
        ylim([0, 70]); 
        
        x_tick_vals = linspace(-18, 18, 7); 
        y_tick_vals = linspace(0, 70, 8);   
        
        xticks(x_tick_vals);
        yticks(y_tick_vals);
        
        x_labels = repmat({''}, size(x_tick_vals));
        x_labels{1} = '-18'; 
        x_labels{end} = '18';
        xticklabels(x_labels);
        
        y_labels = repmat({''}, size(y_tick_vals));
        y_labels{1} = '0'; 
        y_labels{end} = '70';
        yticklabels(y_labels);
        % ---------------------------
        
        hold off;
        
        ax1 = gca;
        pos1 = get(ax1, 'Position');
        
        ax2 = axes('Position', pos1, 'Color', 'none', 'XAxisLocation', 'top', ...
                   'YAxisLocation', 'right', 'HandleVisibility', 'off', 'Visible', 'off');
        
        colormap(ax2, cmap);
        hcb2 = colorbar(ax2, 'Position', [pos1(1) + pos1(3) + 0.05, pos1(2), 0.03, pos1(4)]);
        caxis(ax2, [0, trans_max]);
        
        set(hcb2, 'Ticks', [0, trans_max], 'TickLabels', {'0.0', '1.0'}, ...
                  'FontSize', 14, 'Color', 'k');
        
        ylabel(hcb2, 'Transmission Amplitude (Normalized)', 'FontSize', 16, ...
               'FontWeight', 'bold', 'Color', 'k');
        
        colormap(ax1, gray);
        
        set(ax1, 'Position', [pos1(1), pos1(2), pos1(3) - 0.05, pos1(4)]);
        
        % Customizing the Data Tips
        dcm = datacursormode(gcf);
        set(dcm, 'UpdateFcn', @(obj, event_obj) customDataTip(event_obj, voltage, x, y));

        % fprintf('Transmission mode: %g\n', maxValuesRaw);
        
    %% Produce resolution plots
    [lmaxY, lmaxX] = find(lmax);       % row, col for each local max
    validProm = prom(lmax);            % the prominence values at those same pixels
    
    if isempty(validProm)
        warning('No local maxima found, cannot pick highest-prominence maximum.');
        return;
    end
    
    [~, idxMaxProm] = max(validProm);  % find the one highest-prominence maximum
    
    bestY = lmaxY(idxMaxProm);  
    bestX = lmaxX(idxMaxProm);
    
    % Extract the intensity slices at the chosen maximum location:
    rowIntensity = intensity_rayleigh(bestY, :);
    colIntensity = intensity_rayleigh(:, bestX);

    % colIntensity = intensity_rayleigh(:, int32(size(intensity_rayleigh, 2)./2)); %% experimental line
    
    % Normalize the intensity values to the [0, 1] range:
    normRowIntensity = (rowIntensity - min(rowIntensity)) / (max(rowIntensity) - min(rowIntensity));
    normColIntensity = (colIntensity - min(colIntensity)) / (max(colIntensity) - min(colIntensity));
    
    % Identify local maxima from the normalized signals with their prominences:
    [TF, P] = islocalmax(normRowIntensity, 'FlatSelection', 'first');  % for horizontal slice
    [TFC, PC] = islocalmax(normColIntensity, 'FlatSelection', 'first');  % for vertical slice
    
    % Keep only the few highest-prominence local maxima for the horizontal slice:
    rowMaxIndices = find(TF);
    if numel(rowMaxIndices) > 2
        [~, sortOrder] = sort(P(rowMaxIndices), 'descend');
        topRowIndices = rowMaxIndices(sortOrder(1:2));
        TF_new = false(size(TF));
        TF_new(topRowIndices) = true;
        TF = TF_new;
    end
    
    % Keep only the few highest-prominence local maxima for the vertical slice:
    colMaxIndices = find(TFC);
    if numel(colMaxIndices) > 6
        [~, sortOrder] = sort(PC(colMaxIndices), 'descend');
        topColIndices = colMaxIndices(sortOrder(1:6));
        TFC_new = false(size(TFC));
        TFC_new(topColIndices) = true;
        TFC = TFC_new;
    end

    disp(max(P));
    disp(max(PC));
    
    %% --- Plot the normalized horizontal slice ---
    figure('Color', 'w'); 
    
    plot(x, normRowIntensity, 'b-', x(TF), normRowIntensity(TF), 'r*', 'LineWidth', 2);
    
    xlim([-15, 15]);
    ylim([0, 1]);
    
    xlabel('Distance (cm)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
    ylabel('Normalized Intensity', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
    
    set(gca, 'FontSize', 14, 'FontName', 'Arial', 'LineWidth', 1.5, ...
             'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'TickDir', 'out');
    box on;
    
    x_tick_vals = linspace(-15, 15, 7);
    y_tick_vals = linspace(0, 1, 5);
    
    xticks(x_tick_vals);
    yticks(y_tick_vals);
    
    x_labels = repmat({''}, size(x_tick_vals));
    x_labels{1} = '-15'; 
    x_labels{end} = '15';
    xticklabels(x_labels);
    
    y_labels = repmat({''}, size(y_tick_vals));
    y_labels{1} = '0.0'; 
    y_labels{end} = '1.0';
    yticklabels(y_labels);
    
    %% --- Plot the normalized vertical slice ---
    figure('Color', 'w'); 
    
    plot(y, normColIntensity, 'b-', y(TFC), normColIntensity(TFC), 'r*', 'LineWidth', 2);
    
    xlim([0, 45]);
    ylim([0, 1]);
    
    xlabel('Distance (cm)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
    ylabel('Normalized Intensity', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
    
    set(gca, 'FontSize', 14, 'FontName', 'Arial', 'LineWidth', 1.5, ...
             'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'TickDir', 'out');
    box on;
    
    x_tick_vals_vert = linspace(0, 45, 7);
    y_tick_vals_vert = linspace(0, 1, 5);
    
    xticks(x_tick_vals_vert);
    yticks(y_tick_vals_vert);
    
    x_labels_vert = repmat({''}, size(x_tick_vals_vert));
    x_labels_vert{1} = '0'; 
    x_labels_vert{end} = '45';
    xticklabels(x_labels_vert);
    
    y_labels_vert = repmat({''}, size(y_tick_vals_vert));
    y_labels_vert{1} = '0.0'; 
    y_labels_vert{end} = '1.0';
    yticklabels(y_labels_vert);

    end
    
    function txt = customDataTip(event_obj, voltage, x, y, maxDistance, xResolution, yResolution)
        % Custom data tip function to show the original voltage values and polar coordinates
        pos = get(event_obj, 'Position');
        tolerance = 1e-3; % Define a tolerance for floating-point comparisons
        
        xIndex = find(abs(x - pos(1)) < tolerance);
        yIndex = find(abs(y - pos(2)) < tolerance);
        
        % Retrieve the corresponding voltage value
        voltageValue = voltage(yIndex, xIndex);
        
        % Convert Cartesian to Polar coordinates
        xi = pos(1);
        yi = pos(2);
        
        r = sqrt(xi^2 + yi^2);
        theta = atan2d(yi, xi) - 90; % Adjust for the 90-degree rotation
        
        % Construct the data tip text
        txt = {['X: ', num2str(pos(1))], ...
               ['Y: ', num2str(pos(2))], ...
               ['Voltage: ', num2str(voltageValue)], ...
               ['r: ', num2str(r)], ...
               ['\theta: ', num2str(-theta)]};
    end
end

%% Beamforming %%

function [delayedSig, shiftR1, shiftR2, shiftR3, shiftR4] = apply_delays(ave1, ave2, ave3, ave4, angle, interpNum, interpDen, shiftFactor)
    % Interpolate and downsample to get data with new period
    i_ave1 = interp_downsample(ave1, interpNum, interpDen);
    i_ave2 = interp_downsample(ave2, interpNum, interpDen);
    i_ave3 = interp_downsample(ave3, interpNum, interpDen);
    i_ave4 = interp_downsample(ave4, interpNum, interpDen);

    order = [3, 2, 1, 0];

    % Calculate the required shift for each receiver
    if (angle < 0)
        shiftR4 = circshift(i_ave4, round(order(1) * shiftFactor));
        shiftR3 = circshift(i_ave3, round(order(2) * shiftFactor));
        shiftR2 = circshift(i_ave2, round(order(3) * shiftFactor));
        shiftR1 = circshift(i_ave1, round(order(4) * shiftFactor));
    else
        shiftR1 = circshift(i_ave1, round(order(1) * shiftFactor));
        shiftR2 = circshift(i_ave2, round(order(2) * shiftFactor));
        shiftR3 = circshift(i_ave3, round(order(3) * shiftFactor));
        shiftR4 = circshift(i_ave4, round(order(4) * shiftFactor));
    end

    % Add all the signals together. Take envelope.
    delayedSig = shiftR1 + shiftR2 + shiftR3 + shiftR4;
    
    [delayedSig, ~] = envelope(delayedSig);
end


function newSig = interp_downsample(sig, interpNum, interpDen)
    % Interpolate and downsample
    newSig = downsample(interp(sig, interpNum), interpDen);
end

%% Burst Averaging %%

function [ave1, ave2, ave3, ave4] = average_bursts(rec1, rec2, rec3, rec4, burstSemaphore, samplingRate, dataFileName)
    burstIndices = find(burstSemaphore == 1);
    burstTot = 0;
    ave1Tot = 0;
    ave2Tot = 0;
    ave3Tot = 0;
    ave4Tot = 0;

    burstLength = max(diff(burstIndices));  % Calculate burstLength dynamically.
                                            % 'max' is the correct function
                                            % because it's a dropped bit
                                            % that we are concerned about.

    for i = 1:length(burstIndices)-1
        if burstIndices(i+1) - burstIndices(i) == burstLength
            burstTot = burstTot + 1;
            ave1Tot = ave1Tot + rec1(burstIndices(i):burstIndices(i+1)-1);
            ave2Tot = ave2Tot + rec2(burstIndices(i):burstIndices(i+1)-1);
            ave3Tot = ave3Tot + rec3(burstIndices(i):burstIndices(i+1)-1);
            ave4Tot = ave4Tot + rec4(burstIndices(i):burstIndices(i+1)-1);
        end
    end

    ave1 = ave1Tot / burstTot;
    ave2 = ave2Tot / burstTot;
    ave3 = ave3Tot / burstTot;
    ave4 = ave4Tot / burstTot;
end

%% Data Loading %%

function [rawSignal1, rawSignal2, rawSignal3, rawSignal4, burstSemaphore, angle] = get_receiver_data(dataFileName, samplingRate)
    % Load the data.
    loadedData = load(dataFileName);
    zeroAngle = 90; % change this if needed.
    
    % Grab data and place into variables.
    rawSignal1 = double(loadedData.dataCaptureOut.Rec1);
    rawSignal2 = double(loadedData.dataCaptureOut.Rec2);
    rawSignal3 = double(loadedData.dataCaptureOut.Rec3);
    rawSignal4 = double(loadedData.dataCaptureOut.Rec4);
    angle = mean(double(loadedData.dataCaptureOut.angle)) - zeroAngle; % Adjust angle immediately after loading
    burstSemaphore = loadedData.dataCaptureOut.burst_semaphore;

    % Voltage conversion to millivolts.
    rawSignal1 = (rawSignal1 / (2^12 - 1)) * 1000;
    rawSignal2 = (rawSignal2 / (2^12 - 1)) * 1000;
    rawSignal3 = (rawSignal3 / (2^12 - 1)) * 1000;
    rawSignal4 = (rawSignal4 / (2^12 - 1)) * 1000;
end

%% Transmission Mode Helper Functions %%

function filesTx = sort_files_within_control(filesTx)
    % Get control sequences (fields) from filesTx
    controlSequences = fieldnames(filesTx);
    
    for i = 1:length(controlSequences)
        ctrl = controlSequences{i};
        
        % Extract files for the current control folder
        files = filesTx.(ctrl);
        
        % Initialize an array to hold angles
        angles = zeros(size(files));
        
        % Extract angle from each filename
        for j = 1:length(files)
            filename = files(j).name;
            angleStr = regexp(filename, '-?\d+', 'match');
            angles(j) = str2double(angleStr{end});  % Assumes the angle is the last number in the filename
        end
        
        % Sort files by angle
        [~, sortedIndices] = sort(angles);
        filesTx.(ctrl) = files(sortedIndices);  % Reorder the files in place in filesTx
    end
end

% Helper function to calculate the true maximum for each signal
function trueMax = calculate_segmented_max(signal, b, a, samplingRate, burstSem)
    % Define parameters
    segmentDuration = 10;  % Segment duration in milliseconds
    segmentLength = round(segmentDuration * samplingRate / 1000);  % Segment length in samples

    % Filter the signal
    filteredSignal = filtfilt(b, a, double(signal));

    % Initialize array to store max values for each segment
    maxValues = [];

    % Loop through each 10 ms segment
    for startIdx = 1:segmentLength:length(filteredSignal)-segmentLength+1
        endIdx = startIdx + segmentLength - 1;
        maxValues = [maxValues; max(filteredSignal(startIdx:endIdx))];
    end

    % Calculate the average of the segment-wise max values
    trueMax = mean(maxValues);
end