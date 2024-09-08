%% Clear workspace
% 2024-09-04: Created by Monica Petersen. Note, the pipeline uses the fidall, SPM12, Image Processing toolbox and Nifti toolbox.
clc, clear, close all
%% Load data
% Load waveforms
wfn_b1 = 'C:/Users/hmpet/OneDrive/Skrivebord/MR-centeret/alt-postpros/wavefiles/floret_23Na_fov200_mtx35_arms0p7_hub3_1_intlv120_kdt12_gmax30_smax119_dur15p4.wav'; % Change to you own path
wfn = 'C:/Users/hmpet/OneDrive/Skrivebord/MR-centeret/alt-postpros/wavefiles/cart3D_23Na_fov150_mtx80_nexc6400_kdt84_gmax32_smax119_dur15p5_bal.wav';               % Change to you own path
% Pathname and filename
scriptdir = pwd;
datadir = uigetdir; datadir = [datadir '\'];
cd(datadir)

% Load B1 transmit map
[filename,pathname,~] = uigetfile('*.h5','Load B1 transmit map'); 
fname_b1 = [pathname filename];

% Load 3D Cart sodium images
[filename,pathname,~] = uigetfile('*.h5','Load 3D cart sodium images'); 
fname = [pathname filename];
cd(scriptdir)

%% Recon for 23Na and B1 map
scaling_zero = 1;                                                         % Scaling factor of zero fill for this pipeline dont change it
waveform_info = load_waveform(wfn);
zerofill = scaling_zero*waveform_info.mtx(1);                             % Recon image matrix size

% B1 map
[~,header] = read_p(fname_b1);
lbt_b1 = 0;                                                               % Gaussian linebroadening [Hz] for further prosses dont change
[bb_b1,bbabs_b1] = recon_grid3d(fname_b1,header,wfn_b1,zerofill*[1 1 1],[],lbt_b1);
b = zeros(size(bb_b1,1),size(bb_b1,2),size(bb_b1,4),1,size(bb_b1,3),size(bb_b1,5)); % Setting up the B1 transmit map for BLOSI on 23Na dim(mtx,mtx,#exc,#metabolites,#slice,#coils)

% Seperate the two pulses
if size(bb_b1,5) > 1 % check number of coil elements
    b(:,:,1,:,:,:) = bb_b1(:,:,:,1,:);
    b(:,:,2,:,:,:) = bb_b1(:,:,:,2,:);
else
    b(:,:,1,:,:,:) = bb_b1(:,:,:,1);
    b(:,:,2,:,:,:) = bb_b1(:,:,:,2);
end

[b1map,bbabs_b1map] = blosi_b1map_low_snr(b,header);

% Sodium images
[~,header] = read_p(fname);
lbt = 0;                                                                  % Gaussian linebroadening [Hz] for further prosses dont change
[bb_org,bbabs] = recon_cart(fname,header,wfn,zerofill*[1 1 1],[],lbt);

%% Step 2: denoising

% Define the kernel size
kernel_size = [5 5 5];
  
data_dn = denoise_recursive_tensor(bb_org, kernel_size, 'indices', {1:3 4});  % Denoising function 

% 
bb= mean(data_dn,4);                                                     % Averaging the images into 1

% B1 field correct
b1mapping_adj = fiex3d(b1map);
b1mapping_adj(b1mapping_adj==0) = nan;                                   % avoid dividing with zero
b1mapping_norm = b1mapping_adj./max(b1mapping_adj(:));                   % normalize
bb_b1_corr = bb./b1mapping_norm;                                         % adjust apodised image with b1 correction
bb_b1_corr(isinf(bb_b1_corr)|isnan(bb_b1_corr)) = max(b1mapping_adj(:)); % nan is set to max

% Adjust NaNs in bb_apo_b1_corr
bb_b1_corr(isnan(b1mapping_adj)) = nan;
bb_b1_corr(isnan(b1mapping_adj)) = 0;

bbabs_b1_corr=abs(bb_b1_corr);

% PVC

 disp('PVC start')
 
 iter = 3;
 lb = [10 10];
 bb1_coilC = zeros(size(abs(bb_b1_corr)));
 % bb2_coilC = zeros(size(bbabs2));
 for i = 1:size(bb_b1_corr,3)
     for j = 1:size(bb_b1_corr,4)
         for k = 1:size(bb_b1_corr,5)
             [bb1_coilC(:,:,i,j,k),~]=PVC_sodium_MP(squeeze(bbabs_b1_corr(:,:,i,j,k)),iter,wfn,lb);
         end
     end
 end
 % bbabs= bb1_coilC;
 disp('PVC slut')

 bb1_coilC_single=single(bb1_coilC);
 
% Step 3: Zero fill

 target_size = [160, 160, 160];                                           % doubles the matrix
 upscaled_data = zerofilling(bb1_coilC_single, target_size);

% Visualize post-pros steps
figure ('Name', 'Recon data (bb\_org)'), imagesc_ind3d(abs(mean(bb_org,4))), colormap default,title('Mean of Original Dataset (bb\_org)');
figure ('Name', 'Denoised and B1-Corrected Data (bb\_b1\_corr)'), imagesc_ind3d(abs(bb_b1_corr)),colormap default, title('B1-Corrected Data (bb\_b1\_corr)');
figure ('Name', 'Partial Volume Correction Data (bb1\_coilC\_single)'), imagesc_ind3d(abs(bb1_coilC_single)),colormap default, title('Single Coil-Combined Data (bb1\_coilC\_single)');
figure ('Name', 'Zero filled Data (upscaled\_data)'), imagesc_ind3d(abs(upscaled_data)),colormap default, title('Upscaled Data (upscaled\_data)');

%% Create NifTI

% Step 1: Prepare data (separate real and imaginary components)
real_data = real(abs(upscaled_data));                     % Real part of the data
imag_data = imag(upscaled_data);                          % Imaginary part of the data

% Step 2: Create NIfTI structures for real and imaginary parts
nii_real = make_nii(real_data);                           % Create NIfTI object for real data
nii_imag = make_nii(imag_data);                           % Create NIfTI object for imaginary data

% Step 3: Save the NIfTI files (for both real and imaginary components)
save_nii(nii_real, 'real_part.nii');                      % Save the real part as 'real_part.nii'
save_nii(nii_imag, 'imaginary_part.nii');                 % Save the imaginary part as 'imaginary_part.nii'

% Optionally, you can combine both real and imaginary parts into a single 4D NIfTI
combined_data = cat(4, real_data, imag_data);             % Combine real and imaginary parts
nii_combined = make_nii(combined_data);                   % Create a 4D NIfTI
save_nii(nii_combined, '23Na_postpros_mtx160x160.nii');   % Save the 4D file


%% SPM12 

% t1 seq
spm('defaults', 'fmri');
spm_jobman('initcfg');

matlabbatch = {};

% Step 1: Prompt the user to select the T1-weighted NIfTI file
[t1_file, t1_path] = uigetfile('*.nii', 'Select the T1-weighted NIfTI file');
if isequal(t1_file, 0)
    disp('User canceled the file selection.');
    return;
end
t1_full_path = fullfile(t1_path, t1_file);

% Step 2: Prompt the user to select the folder where tissue classes and deformation fields are saved
output_folder = uigetdir(pwd, 'Select the folder to save tissue classes and deformation fields');
if isequal(output_folder, 0)
    disp('User canceled the folder selection.');
    return;
end

% Bias correction settings
matlabbatch{1}.spm.spatial.preproc.channel.vols = {t1_full_path};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];                 % Save bias field and bias-corrected image

% Specify tissue probability maps (TPMs) and settings
tpm_path = 'C:/MatLab/SPM12/spm12/spm12/tpm/TPM.nii';                     % Change to the correct folder on your PC
for i = 1:6
    matlabbatch{1}.spm.spatial.preproc.tissue(i).tpm = {[tpm_path ',' num2str(i)]};
    matlabbatch{1}.spm.spatial.preproc.tissue(i).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(i).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(i).warped = [0 0];          % Not requesting warped tissue images here
end

% Warping settings
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = 4;
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';                   % can be changed to 'ecc' or 'mi'
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;

% Save the deformation field (inverse and forward)
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1]; % Save deformation fields

% Paths for tissue classes (c1, c2, c3) and forward deformation field
tissues = {
    fullfile(output_folder, ['c1_' t1_file]), 
    fullfile(output_folder, ['c2_' t1_file]), 
    fullfile(output_folder, ['c3_' t1_file])
};

deformation_field = fullfile(output_folder, ['y_' t1_file]);

% Run the batch
spm_jobman('run', matlabbatch);


%% Coregistrate t1 to atlas 

spm('defaults', 'fmri');
spm_jobman('initcfg');

% Choose path to atlas
[atlas_file, atlas_path] = uigetfile('*.nii', 'Select atlas NIfTI file');
if isequal(atlas_file, 0)
    disp('User canceled the file selection.');
    return;
end
atlas_full_path = fullfile(atlas_path, atlas_file);

% Coregistration batch
matlabbatch = {};

% Coregister Estimate & Reslice
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {sprintf('%s,1', t1_full_path)};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {sprintf('%s,1', atlas_full_path)};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'ecc';      % otherwise use 'nmi' or 'mi'
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

% Run the job
spm_jobman('run', matlabbatch);

%% Coregister and segment 23Na images

spm('defaults', 'fmri');
spm_jobman('initcfg');

% Coregister the 23Na images to the T1-weighted image, load image 23Na
[na_file, na_path] = uigetfile('*.nii', 'Select the post-pros (x-) 23Na .nii');
if isequal(na_file, 0)
    disp('User canceled the file selection.');
    return;
end
na_full_path = fullfile(na_path, na_file);

[t1_file, t1_path] = uigetfile('*.nii', 'Select the (x-) t1 .nii');
if isequal(t1_file, 0)
    disp('User canceled the file selection.');
    return;
end
t1_full_path = fullfile(t1_path, t1_file);

% Coregistration batch
matlabbatch = {};

% Coregister Estimate & Reslice
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {sprintf('%s,1', t1_full_path)};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {sprintf('%s,1',na_full_path )};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'ecc';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [2 1];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

% Run the job
spm_jobman('run', matlabbatch);

%% OBS Inspect the coregistration manually in ITK-snap
% Sometimes the alignment is not quiet satisfyring even after changing the
% parameters within the coreg therefore check it visualy and correct the
% position if needed. 

% A guide to correction is found in labbook under Monica -->

%% Apply mask to 23Na and measure concentration in WM, GM, and CSF

% Load necessary SPM functions
spm('defaults', 'FMRI');

% Paths to your coregistered 23Na image and masks
% Prompt user to select the Sodium (23Na) image
[Na_img_file, Na_img_path] = uigetfile('*.nii', 'Select the coregistred (rx) Sodium (23Na) NIfTI Image');
if isequal(Na_img_file, 0)
    disp('User canceled the selection of the Sodium image.');
    return;
else
    Na_img_full_path = fullfile(Na_img_path, Na_img_file);
end

% Prompt user to select the Gray Matter (GM) mask
[GM_mask_file, GM_mask_path] = uigetfile('*.nii', 'Select the c1- Gray Matter (GM) Mask');
if isequal(GM_mask_file, 0)
    disp('User canceled the selection of the GM mask.');
    return;
else
    GM_mask_full_path = fullfile(GM_mask_path, GM_mask_file);
end

% Prompt user to select the White Matter (WM) mask
[WM_mask_file, WM_mask_path] = uigetfile('*.nii', 'Select the c2- White Matter (WM) Mask');
if isequal(WM_mask_file, 0)
    disp('User canceled the selection of the WM mask.');
    return;
else
    WM_mask_full_path = fullfile(WM_mask_path, WM_mask_file);
end

% Prompt user to select the Cerebrospinal Fluid (CSF) mask
[CSF_mask_file, CSF_mask_path] = uigetfile('*.nii', 'Select the c3- Cerebrospinal Fluid (CSF) Mask');
if isequal(CSF_mask_file, 0)
    disp('User canceled the selection of the CSF mask.');
    return;
else
    CSF_mask_full_path = fullfile(CSF_mask_path, CSF_mask_file);
end

% Load the images
Na_img = spm_vol(Na_img_full_path);
Na_data = spm_read_vols(Na_img);
GM_mask = spm_vol(GM_mask_full_path);
GM_data = spm_read_vols(GM_mask);
WM_mask = spm_vol(WM_mask_full_path);
WM_data = spm_read_vols(WM_mask);
CSF_mask = spm_vol(CSF_mask_full_path);
CSF_data = spm_read_vols(CSF_mask);

% Threshold the masks to create binary masks (you can adjust the threshold)
GM_bin = GM_data > 0.5;
WM_bin = WM_data > 0.5;
CSF_bin = CSF_data > 0.5;

% Calculate mean values for Na within each mask
GM_Na_values = Na_data(GM_bin);
WM_Na_values = Na_data(WM_bin);
CSF_Na_values = Na_data(CSF_bin);

% Compute mean sodium values for each ROI (GM, WM, CSF)
mean_GM_Na = mean(GM_Na_values, 'omitnan');
mean_WM_Na = mean(WM_Na_values, 'omitnan');
mean_CSF_Na = mean(CSF_Na_values, 'omitnan');

% Display the results
fprintf('Mean Sodium in GM: %f\n', mean_GM_Na);
fprintf('Mean Sodium in WM: %f\n', mean_WM_Na);
fprintf('Mean Sodium in CSF: %f\n', mean_CSF_Na);

% Measure 23Na concentration in specific ROIs (excluding CSF)

% Load the NeuroMorphometrics Atlas NIfTI file
[atlas_file, atlas_path] = uigetfile('*.nii', 'Select the coregistred atlas (r- ) labels_Neuromorphometrics.nii');
if isequal(atlas_file, 0)
    disp('User canceled the selection of the atlas.');
    return;
else
    atlas_full_path = fullfile(atlas_path, atlas_file);
end

atlasData = niftiread(atlas_full_path);

% List of ROI indices for the desired brain regions (CSF removed)
roiIndices = [4, 11, 31, 32, 47, 48, 51, 52, 59, 60]; 

% ROI Names for Reference (without CSF)
roiNames = { ...
    '3rd Ventricle', ...
    '4th Ventricle', ...
    'Right Amygdala', ...
    'Left Amygdala', ...
    'Right Hippocampus', ...
    'Left Hippocampus', ...
    'Right Lateral Ventricle', ...
    'Left Lateral Ventricle', ...
    'Right Thalamus Proper', ...
    'Left Thalamus Proper'};

% Initialize results
concentrationResults = zeros(length(roiIndices), 1);

% Loop through each ROI and calculate the mean 23Na concentration
for i = 1:min(length(roiIndices), length(roiNames))
    roiMask = (atlasData == roiIndices(i));
    roiConcentration = Na_data(roiMask);
    meanConcentration = mean(roiConcentration, 'omitnan');  % Calculate mean concentration excluding NaN
    concentrationResults(i) = meanConcentration;

    fprintf('Mean 23Na Concentration in %s: %.2f\n', roiNames{i}, meanConcentration);
end

%% Add WM, GM, and CSF values to the ROI results and names

roiNamesWithWMGMCSF = [{'Gray Matter', 'White Matter', 'Cerebrospinal Fluid'}, roiNames];
concentrationResultsWithWMGMCSF = [mean_GM_Na; mean_WM_Na; mean_CSF_Na; concentrationResults];

% Ensure both variables have the same length
if length(roiNamesWithWMGMCSF) == length(concentrationResultsWithWMGMCSF)
    % Create a table with region names and sodium concentrations
    resultsTable = table(roiNamesWithWMGMCSF', concentrationResultsWithWMGMCSF, ...
                         'VariableNames', {'Region', 'Mean23NaConcentration'});

    % Save the results to a CSV file
    writetable(resultsTable, '23Na_signal_intensity_Results_WM_GM_CSF_ROI.csv');
    disp('23Na signal intesity calculation completed and saved to 23Na_SignalIntensity_Results_WM_GM_CSF_ROI.csv');
else
    disp('Error: Mismatch in the number of regions and concentrations.');
end

%% For quantification of Na concentration within mmol/L go to ITK-snap and manual segnmentate phantoms 

% Use label 1 (red) for segmentate 0.4g NaCl phantom
% Use label 2 (green) for segmentate 0.1 NaCl phantom

%% Calculation of Na concentration during calibration curve

[Seg_file, Seg_path] = uigetfile('*.nii', 'Select your phantom segmentation');
if isequal(Na_img_file, 0)
    disp('User canceled the selection of the Sodium image.');
    return;
else
    Seg_full_path = fullfile(Na_img_path, Na_img_file);
end

% Calculation of Na concentration during calibration curve

% Load the Sodium image and segmentation data
Na_info = niftiinfo(Na_img_full_path);
Na_data = niftiread(Na_img_full_path);

Seg_info = niftiinfo(Seg_full_path);
Seg_data = niftiread(Seg_full_path);

% Ensure data is in double precision
Na_data = double(Na_data);
Seg_data = double(Seg_data);

% Extract the unique labels from the segmentation
labels = unique(Seg_data);
labels(labels == 0) = [];  % Remove the background label (0), if present

% Display found labels for debugging
disp('Unique labels found in the segmentation:');
disp(labels);

% Check if the number of labels is correct (should be exactly 2)
if numel(labels) ~= 2
    error('Unexpected number of labels. There should be exactly two labels in the segmentation.');
end

% Assign labels (assumed 1 for high concentration and 2 for low concentration)
phantom4g_label = 1;
phantom1g_label = 2;

% Sodium concentrations for each phantom (from calibration data)
na_concentration_low = 34.22138;   % mmol/L for 01phantom (low concentration)
na_concentration_high = 136.8855;  % mmol/L for 04phantom (high concentration)

% Extract mean signal intensity for each label
phantom1_signal = Na_data(Seg_data == phantom1g_label);
mean_phantom1_signal = mean(phantom1_signal);

phantom2_signal = Na_data(Seg_data == phantom4g_label);
mean_phantom2_signal = mean(phantom2_signal);

% Define the phantom data points
phantom_concentrations = [na_concentration_low; na_concentration_high];   % x-values (mmol/L)
mean_signals = [mean_phantom1_signal; mean_phantom2_signal];              % y-values (signal intensity)

% Fit the linear model: y = ax + b
p_linear = polyfit(phantom_concentrations, mean_signals, 1);

% Display the fitted coefficients
a = p_linear(1);  % Slope
b = p_linear(2);  % Intercept
fprintf('Fitted model: y = %.4f * x + %.4f\n', a, b);

% Generate a range of sodium concentrations for plotting
concentration_range = linspace(min(phantom_concentrations) * 0.8, max(phantom_concentrations) * 1.2, 100);

% Calculate corresponding signal intensities using the fitted model
signal_intensity_linear = polyval(p_linear, concentration_range);

% Plot the calibration curve
figure;
plot(concentration_range, signal_intensity_linear, '-r', 'LineWidth', 2); % Linear fit
hold on;

% Plot the actual phantom data points
plot(phantom_concentrations, mean_signals, 'ok', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

% Add labels and legend
xlabel('Na Concentration (mmol/L)');
ylabel('Mean Signal Intensity');
title('Calibration Curve for Sodium Concentration');
legend('Fitted Linear Model', 'Phantom Data', 'Location', 'best');
grid on;
hold off;

%%
% Prompt user to select the CSV file containing the mean signal intensity data
[mean_signal_file, mean_signal_path] = uigetfile('*.csv', 'Select the 23Na Concentration Results CSV file');
if isequal(mean_signal_file, 0)
    disp('User canceled the file selection.');
    return;
else
    mean_signal_full_path = fullfile(mean_signal_path, mean_signal_file);
end

% Load the mean signal intensity data
mean_signal_data = readtable(mean_signal_full_path); 

% Extract the mean signal intensities for each ROI
mean_signal_intensities = mean_signal_data.Mean23NaConcentration; % Assuming the column is named 'MeanSignalIntensity'

% Define ROIs
ROIs = {'GM',...
        'WM',...
        'CSF',...
        '3rd Ventricle', ...
        '4th Ventricle', ...
        'Right Amygdala', ...
        'Left Amygdala', ...
        'Right Hippocampus', ...
        'Left Hippocampus', ...
        'Right Lateral Ventricle', ...
        'Left Lateral Ventricle', ...
        'Right Thalamus Proper', ...
        'Left Thalamus Proper'};

% Initialize an array to store the sodium concentrations
Na_concentrations = zeros(length(ROIs), 1);

% Display the initialized Na_concentrations array
disp(Na_concentrations);


% Initialize an array to store the sodium concentrations
Na_concentrations = zeros(length(ROIs), 1);

% Calculate sodium concentration for each ROI
for i = 1:length(ROIs)
    mean_signal = mean_signal_intensities(i); % Get the mean signal intensity for the ROI
    Na_concentrations(i) = (mean_signal - b) / a; % Apply the calibration curve equation
end

% Display the sodium concentrations for each ROI
for i = 1:length(ROIs)
    fprintf('%s Na concentration: %.2f mmol/L\n', ROIs{i}, Na_concentrations(i));
end

% Create a table with the ROIs and their sodium concentrations
T = table(ROIs', Na_concentrations, 'VariableNames', {'ROI', 'Na_Concentration_mmol_L'});

% Save the table to a CSV file
writetable(T, 'sodium_concentrations.csv');

% Confirm the file has been saved
disp('Data has been saved to sodium_concentrations.csv');