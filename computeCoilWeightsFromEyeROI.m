function [weightsNorm, eyeMask2D, irMask2D, coilOrder, info] = computeCoilWeightsFromEyeROI(xrmsOrPath, C, varargin)
% computeCoilWeightsFromEyeROI
% Rank/weight coils by their signal magnitude inside a manually drawn eye ROI.
%
% xrmsOrPath can be:
%   (1) a path to a .mat file containing a 3D (or 2D) reference volume, OR
%   (2) the reference volume itself (xrms), size ideally [Nx Ny Nz]
%
% The reference volume will be resampled to match C spatial size [Nx Ny Nz]
% (undersampled if larger; generally resized either way for safety).
%
% INPUTS
%   xrmsOrPath : [] | char/string path | numeric array (2D/3D)
%               If [] or not provided, we try to derive a reference from C (RMS across coils).
%   C          : coil images, size [Nx, Ny, Nz, Ncoil]
%
% OPTIONAL NAME-VALUE PAIRS
%   'Verbose'      : true/false (default: false)
%   'ZRange'       : vector of slices (default: 1:Nz)
%   'NumTop'       : number of top coils to visualize (default: 10)
%   'SliceIndex'   : slice index for ROI drawing (default: mid-slice)
%   'RefMode'      : 'xrms' | 'crms' (default: 'xrms' if provided else 'crms')
%   'Interp'       : 'nearest' | 'linear' | 'cubic' (default: 'linear') for ref resampling
%
% OUTPUTS
%   weightsNorm : normalized coil weights (sum = 1), size [Ncoil,1]
%   eyeMask2D   : drawn ROI mask in 2D (Nx-by-Ny logical)
%   coilOrder   : coil indices sorted by descending weight
%   info        : struct with details (refUsed, zRange, rawWeights, etc.)
% 
%  USAGES:
% [w, mask, order] = computeCoilWeightsFromEyeROI(xrms, C, 'Verbose', true);
% [w, mask, order] = computeCoilWeightsFromEyeROI('/path/xrms.mat', C, 'ZRange', 12:36);
% [w, mask, order] = computeCoilWeightsFromEyeROI([], C);
% -------------------- Defaults / parse inputs --------------------
if nargin < 2
    error('Usage: computeCoilWeightsFromEyeROI(xrmsOrPath, C, ...)');
end

p = inputParser;
p.addParameter('Verbose', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('ZRange', [], @(x)isnumeric(x) && isvector(x));
p.addParameter('NumTop', 10, @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('SliceIndex', [], @(x)isnumeric(x) && isscalar(x));
p.addParameter('RefMode', '', @(s)ischar(s) || isstring(s));
p.addParameter('Interp', 'linear', @(s)ischar(s) || isstring(s));
p.parse(varargin{:});

verbose    = logical(p.Results.Verbose);
zRangeUser = p.Results.ZRange;
numTop     = round(p.Results.NumTop);
sliceIndex = p.Results.SliceIndex;
refMode    = lower(string(p.Results.RefMode));
interpMode = lower(string(p.Results.Interp));

% -------------------- Basic checks --------------------
if ndims(C) ~= 4
    error('C must be 4D: [Nx Ny Nz Ncoil].');
end
[nx, ny, nz, nCoils] = size(C);

% -------------------- Acquire reference volume (xrms) --------------------
% If xrmsOrPath is empty/missing -> derive a reference from C (RMS across coils)
if isempty(xrmsOrPath)
    xrms = deriveRefFromC(C);
    refUsed = "crms";
else
    if ischar(xrmsOrPath) || isstring(xrmsOrPath)
        % load from path; use first variable
        S = load(xrmsOrPath);
        fn = fieldnames(S);
        if isempty(fn); error('MAT file contains no variables.'); end
        xrms = S.(fn{1});
        refUsed = "xrms(path)";
    elseif isnumeric(xrmsOrPath) || islogical(xrmsOrPath)
        xrms = xrmsOrPath;
        refUsed = "xrms(var)";
    else
        error('xrmsOrPath must be [] or a path or a numeric array.');
    end
end

% Decide refMode (display source)
if refMode == ""
    if contains(refUsed, "xrms")
        refMode = "xrms";
    else
        refMode = "crms";
    end
end

% Ensure xrms is 3D
if ndims(xrms) == 2
    xrms = reshape(xrms, size(xrms,1), size(xrms,2), 1);
elseif ndims(xrms) > 3
    error('Reference volume must be 2D or 3D. Got ndims=%d.', ndims(xrms));
end

% -------------------- Resample reference to match C size [Nx Ny Nz] --------------------
% Requirement: xrms spatial size must match C per-coil spatial size.
% If xrms is larger (common), it should be undersampled; in general we resize safely either way.
xrmsRes = resampleVolumeToMatch(xrms, [nx ny nz], interpMode);

% Pick slice for ROI drawing
if isempty(sliceIndex)
    sliceIndex = max(1, round(nz/2));
end
sliceIndex = max(1, min(nz, sliceIndex));

% -------------------- Verbose: show ref --------------------
if verbose
    figure;
    imagesc(abs(xrmsRes(:,:,sliceIndex))); axis image off; colormap gray;
    title(sprintf('Reference for ROI (%s), slice %d', refMode, sliceIndex));
end

% -------------------- Draw ROI on reference slice --------------------
figure;
imagesc(abs(xrmsRes(:,:,sliceIndex))); axis image off; colormap gray;
title('Draw ROI around the eyes, double-click to finish');
eyeMask2D = logical(roipoly);

imagesc(abs(xrmsRes(:,:,sliceIndex))); axis image off; colormap gray;
title('Draw ROI around the interference region, double-click to finish');
irMask2D = logical(roipoly);

if verbose
    figure;
    imagesc(eyeMask2D); axis image off; colormap gray;
    title('Eye ROI mask (2D)');
    figure;
    imagesc(abs(xrmsRes(:,:,sliceIndex)) .* eyeMask2D); axis image off; colormap gray;
    title('ROI applied to reference slice');
    figure;
    imagesc(irMask2D); axis image off; colormap gray;
    title('IR mask (2D)');
    figure;
    imagesc(abs(xrmsRes(:,:,sliceIndex)) .* irMask2D); axis image off; colormap gray;
    title('IR applied to reference region');
end

% -------------------- Decide zRange --------------------
if isempty(zRangeUser)
    zRange = 1:nz;
else
    zRange = unique(round(zRangeUser(:)'));
    zRange = zRange(zRange >= 1 & zRange <= nz);
    if isempty(zRange)
        error('Provided ZRange is empty after clipping to [1..%d].', nz);
    end
end

% -------------------- Compute coil weights in ROI --------------------
mask3D = repmat(eyeMask2D, [1, 1, numel(zRange)]);
maskIdx = mask3D(:);

rawWeights = zeros(nCoils, 1, 'double');
for coil = 1:nCoils
    coilMag = abs(C(:,:,zRange,coil));
    v = coilMag(:);
    rawWeights(coil) = sum(v(maskIdx));
end

% Normalize safely
s = sum(rawWeights);
if s <= 0 || ~isfinite(s)
    warning('Sum of rawWeights is non-positive or invalid. Returning uniform weights.');
    weightsNorm = ones(nCoils,1) / nCoils;
else
    weightsNorm = rawWeights / s;
end

% Sort coils
[~, coilOrder] = sort(weightsNorm, 'descend');

% -------------------- Verbose plots & inspection --------------------
if verbose
    figure;
    bar(weightsNorm(coilOrder), 'LineWidth', 1.2);
    xlabel('Rank (1 = highest)'); ylabel('Normalized weight');
    title('Coil contributions within eye ROI (sorted)');
    grid on;

    disp('Coils sorted by weight (highest to lowest):');
    disp(coilOrder(:)');

    numTop = min(numTop, nCoils);
    topIdx = coilOrder(1:numTop);

    figure;
    bar(weightsNorm, 'FaceColor', [0.7 0.7 0.7]); hold on;
    bar(topIdx, weightsNorm(topIdx), 'FaceColor', [0.9 0.3 0.3]);
    xlabel('Coil index'); ylabel('Normalized weight');
    title(sprintf('Top %d coils (indices: %s)', numTop, num2str(topIdx')));
    legend('All coils', 'Top contributors');
    grid on;

    if exist('bmImage','file') == 2
        for k = 1:numTop
            bmImage(C(:,:,:,topIdx(k)));
            title(sprintf('Coil %d', topIdx(k)));
        end
    end
end

% -------------------- Info struct --------------------
info = struct();
info.refUsed      = refUsed;
info.refMode      = refMode;
info.refSizeIn    = size(xrms);
info.refSizeMatch = size(xrmsRes);
info.zRange       = zRange;
info.rawWeights   = rawWeights;

end

% ===== Helper: derive reference from C =====
function ref = deriveRefFromC(C)
% Robust reference image for drawing ROI: RMS magnitude across coils
Cmag2 = abs(C).^2;
ref = sqrt(mean(Cmag2, 4));
end

% ===== Helper: resample a 2D/3D volume to match target size =====
function Vout = resampleVolumeToMatch(Vin, targetSize, interpMode)
% Resample Vin to targetSize = [Nx Ny Nz]. Uses imresize3 if available.
targetSize = double(targetSize);
Vin = double(Vin);

% Make sure target has 3 elements
if numel(targetSize) ~= 3
    error('targetSize must be [Nx Ny Nz].');
end

inSize = size(Vin);
if numel(inSize) == 2
    inSize(3) = 1;
end

% If already matches, return
if isequal(inSize(1:3), targetSize)
    Vout = Vin;
    return;
end

% Use imresize3 when available (recommended)
if exist('imresize3','file') == 2
    Vout = imresize3(Vin, targetSize, char(interpMode));
else
    % Fallback: resize each slice in-plane + then resize along z by interp1
    % 1) resize x-y for each z
    tmp = zeros(targetSize(1), targetSize(2), inSize(3));
    for k = 1:inSize(3)
        tmp(:,:,k) = imresize(Vin(:,:,k), targetSize(1:2), char(interpMode));
    end
    % 2) resize along z
    if inSize(3) == targetSize(3)
        Vout = tmp;
    else
        zIn  = linspace(1, inSize(3), inSize(3));
        zOut = linspace(1, inSize(3), targetSize(3));
        Vout = zeros(targetSize(1), targetSize(2), targetSize(3));
        for i = 1:targetSize(1)
            for j = 1:targetSize(2)
                Vout(i,j,:) = interp1(zIn, squeeze(tmp(i,j,:)), zOut, char(interpMode), 'extrap');
            end
        end
    end
end

end