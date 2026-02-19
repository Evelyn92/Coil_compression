%% ========================================================================
% Coil Selection & Compression 
% This script cannot run, it's just a tutorial
%  ========================================================================
%  This script demonstrates:
%   1) How to compute ROI-based coil weights
%   2) How to compute ROI / interference-region GEVD weights (ROVir)
%   https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28706
%   3) How to use them inside coilCompressCore
%   4) How to prepare data for Mitosius reconstruction
%
%  Required inputs:
%     xrms  - RMS image
%     x0    - nCh-elements cell, (the input to calculate xrms)
%     C     - original Sensitivity maps
%     reader - raw data reader object
%     CfilePath - path to saved sensitivity maps
% ========================================================================
% PS: I combined moco with coil compression in mitosius
% in recon_mitosius_moco_cc.m, if interested (I guess you all have yours so
% maybe not :) ).
%% ========================================================================
%% 0.0) No need to compute weights for global SVD (case 0 below)
% ========================================================================

%% ========================================================================
%% 0.2) Compute Coil Weights for Eye ROI coil selection (case 2 below)
% ========================================================================
% Computes coil importance weights based on an eye ROI.
% See function comments for required inputs.

[w, mask, irMask2D, order] = computeCoilWeightsFromEyeROI( ...
    xrms, C, ...
    'Verbose', true, ...
    'ZRange', 12:36);

save('w_eye_selection.mat', 'w');
disp('Saved: w_eye_selection.mat');


%% ========================================================================
%% 0.3) Compute ROI / Interference-Region GEVD (case 3 below)
% ========================================================================

% Computes GEVD-based virtual coil weights.
% See function comments for detailed parameter explanation.

[W, lambda, masksOut_rovir] = computeROVirFromEyeROI( ...
    x0, [], [], [], ...
    'IRMode', 'complement_object', ...
    'DilateROIPixels', 2, ...
    'ObjectThresh', 0.05);

save('w_eye_rovir.mat', 'W');
disp('Saved: w_eye_rovir.mat');


%% ========================================================================
%% 1) Load Raw Data and Sensitivity Maps
% ========================================================================

% --- Read raw data (filters: nShotOff + SI enabled)
y_tot = reader.readRawData(true, true);

% --- Load sensitivity maps
S = load(CfilePath, 'C');
C = S.C;
disp(['Sensitivity maps loaded from: ', CfilePath]);


%% ========================================================================
%% 2) Normalize Raw Data (if necessary)
% ========================================================================

% normalize_val can be computed from bmMathilda if needed.
% Here we assume it already exists.

if mean(abs(y_tot(:))) < 1
    y_tot = y_tot / normalize_val;
    disp('y_tot normalized by normalize_val.');
else
    disp('y_tot appears already normalized.');
end


%% ========================================================================
%% 3) Coil Compression Methods
% ========================================================================
% method 0: Global SVD
% method 2: ROI-weighted Top-N selection
% method 3: ROVir / GEVD virtual coils
% ========================================================================

nChTarget = 25;   % Number of compressed channels
saveDiagnostics = 1;   % Save verbose plots
mDir = someFolder;     % Output directory for diagnostics


%% ------------------------------------------------------------------------
%% 3.0) Case 0 — Global SVD
% ------------------------------------------------------------------------

[y_tot_comp, C_comp, info] = coilCompressCore( ...
    y_tot, C, ...
    'method', 0, ...
    'nChCompressed', nChTarget, ...
    'saveDiagnostics', saveDiagnostics, ...
    'mDir', mDir);


%% ------------------------------------------------------------------------
%% 3.2) Case 2 — ROI-weighted Top-N Selection
% ------------------------------------------------------------------------

[y_tot_comp, C_comp, info] = coilCompressCore( ...
    y_tot, C, ...
    'method', 2, ...
    'nChCompressed', nChTarget, ...
    'weightsNormPath', 'path/to/w_eye_selection.mat', ...
    'saveDiagnostics', saveDiagnostics, ...
    'mDir', mDir);


%% ------------------------------------------------------------------------
%% 3.3) Case 3 — ROVir / GEVD Virtual Coils
% ------------------------------------------------------------------------

[y_tot_comp, C_comp, info] = coilCompressCore( ...
    y_tot, C, ...
    'method', 3, ...
    'nChCompressed', nChTarget, ...
    'weightsNormPath', 'path/to/w_eye_rovir.mat', ...
    'saveDiagnostics', saveDiagnostics, ...
    'mDir', mDir);


%% ========================================================================
%% 4) Save Compressed Sensitivity Maps
% ========================================================================

save("C_comp.mat", 'C_comp');
disp('Saved: C_comp.mat');


%% ========================================================================
%% 5) Prepare Data for Mitosius Reconstruction
% ========================================================================

% --- Split compressed raw data
[y, t] = bmMitosis(y_tot_comp, t_tot, tempMask);

% --- Reformat for reconstruction
y  = bmPermuteToCol(y);
ve = bmVolumeElement(t, 'voronoi_full_radial3');

% --- Create Mitosius files
bmMitosius_create(mDir, y, t, ve);

disp('Mitosius files saved to:');
disp(mDir);


%% ========================================================================
%% 6) Reconstruction
% ========================================================================
% IMPORTANT:
% Replace original C with C_comp in your reconstruction pipeline.
%
% Example:
%     recon = bmSteva(y, C_comp, ...);
%
% ========================================================================