function [W, lambda, masksOut] = computeROVirFromEyeROI(x0, mask, irMask2D, xrms, varargin)
debugFlag = 0;
% computeROVirFromEyeROI
%   Compute ROVir (GEVD) virtual coil weights from an eye ROI and an interference region.
%
% Options (name-value):
%   'IRMode'          : 'manual' (default) | 'complement' | 'complement_object'
%   'DilateROIPixels' : nonnegative integer (default 0)
%   'ObjectThresh'    : in [0,1], relative threshold on abs(xrms) for object mask (default 0.05)
%   'EpsFactor'       : regularization factor multiplier (default 1e-6)
% Usage:
% [W, lambda] = computeROVirFromEyeROI(x0, mask, irMask2D);

% [W, lambda] = computeROVirFromEyeROI(x0, mask, [], xrms, ...
    % 'IRMode','complement', 'DilateROIPixels', 6);

% [W, lambda] = computeROVirFromEyeROI(x0, mask, [], xrms, ...
    % 'IRMode','complement_object', 'DilateROIPixels', 6, 'ObjectThresh', 0.05);
% ---------------- Options ----------------
p = inputParser;
p.addParameter('IRMode', 'manual', @(s) ischar(s) || isstring(s));
p.addParameter('DilateROIPixels', 2, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('ObjectThresh', 0.05, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('EpsFactor', 1e-6, @(x) isnumeric(x) && isscalar(x) && x>0);
p.parse(varargin{:});
opt = p.Results;

opt.IRMode = lower(string(opt.IRMode));

% ---------------- Load x0 if it is a file path ----------------
if ischar(x0) || isstring(x0)
    x0path = char(x0);
    if ~exist(x0path, 'file')
        error('computeROVirFromEyeROI:FileNotFound', 'x0 path not found: %s', x0path);
    end
    s = load(x0path);
    if ~isfield(s, 'x0')
        error('computeROVirFromEyeROI:MissingVar', 'MAT file does not contain variable ''x0'': %s', x0path);
    end
    x0 = s.x0;
end

% ---------------- Convert x0 -> g_cal 4D ----------------
g_cal = x0;
if iscell(g_cal)
    try
        g_cal = cat(4, g_cal{:});   % [Nx,Ny,Nz,nCh]
    catch ME
        error('computeROVirFromEyeROI:BadCell', ...
            'Failed to cat(4, x0{:}). Ensure each cell contains same-size [Nx,Ny,Nz].\n%s', ME.message);
    end
end
if ndims(g_cal) ~= 4
    error('computeROVirFromEyeROI:BadX0', 'x0 must be cell of 3D or 4D [Nx,Ny,Nz,nCh].');
end

[Nx, Ny, Nz, nCh] = size(g_cal);

% ---------------- Get a 2D view of xrms for interactive drawing / object mask ----------------
xrms2D = [];
sliceIndex = 1;
if ~isempty(xrms)
    if ndims(xrms) == 3
        sliceIndex = round(size(xrms,3)/2);
        xrms2D = xrms(:,:,sliceIndex);
    else
        xrms2D = xrms;
    end
end

% ---------------- If ROI missing: interactive draw (needs xrms) ----------------
if isempty(mask)
    if isempty(xrms2D)
        error(['mask is empty but xrms is not provided. ', ...
               'Tip: run check_orient_xrms to generate xrms, then call this function again.']);
    end
    figure; imagesc(abs(xrms2D)); axis image off; colormap gray;
    title('Draw ROI around the eyes, double-click to finish');
    mask = logical(roipoly);
    close(gcf);
    warning('computeROVirFromEyeROI:InteractiveROI', ...
        'ROI mask was empty. Used interactive roipoly (sliceIndex=%d).', sliceIndex);
end

% ---------------- Validate ROI mask size ----------------
if ~isequal(size(mask), [Nx, Ny])
    error('computeROVirFromEyeROI:MaskSizeMismatch', ...
        'ROI mask size must match [Nx,Ny]=[%d,%d], but got [%d,%d].', ...
        Nx, Ny, size(mask,1), size(mask,2));
end

% ---------------- Optionally dilate ROI (2D) ----------------
roi2D = logical(mask);
if opt.DilateROIPixels > 0
    se = strel('disk', opt.DilateROIPixels, 0);
    roi2D = imdilate(roi2D, se);
end

% ---------------- Build IR mask according to IRMode ----------------
switch opt.IRMode
    case "manual"
        % If IR missing: interactive draw (needs xrms)
        if isempty(irMask2D)
            if isempty(xrms2D)
                error(['irMask2D is empty but xrms is not provided. ', ...
                       'Tip: run check_orient_xrms to generate xrms, then call this function again.']);
            end
            figure; imagesc(abs(xrms2D)); axis image off; colormap gray;
            title('Draw ROI around the interference region, double-click to finish');
            irMask2D = logical(roipoly);
            close(gcf);
            warning('computeROVirFromEyeROI:InteractiveIR', ...
                'IR mask was empty. Used interactive roipoly (sliceIndex=%d).', sliceIndex);
        end

        if ~isequal(size(irMask2D), [Nx, Ny])
            error('computeROVirFromEyeROI:IRMaskSizeMismatch', ...
                'IR mask size must match [Nx,Ny]=[%d,%d], but got [%d,%d].', ...
                Nx, Ny, size(irMask2D,1), size(irMask2D,2));
        end
        ir2D = logical(irMask2D);

    case "complement"
        % IR = everything except (dilated) ROI
        ir2D = ~roi2D;

    case "complement_object"
        % IR = object mask (from xrms) minus (dilated) ROI
        if isempty(xrms2D)
            error(['IRMode=complement_object requires xrms. ', ...
                   'Tip: run check_orient_xrms to generate xrms, then call this function again.']);
        end
        mag = abs(xrms2D);
        thr = opt.ObjectThresh * max(mag(:));
        obj2D = mag > thr;
        % optional cleanup
        obj2D = imfill(obj2D, 'holes');
        obj2D = bwareaopen(obj2D, 50);
        ir2D = obj2D & ~roi2D;

    otherwise
        error('computeROVirFromEyeROI:BadIRMode', ...
            'Unknown IRMode: %s. Use ''manual'', ''complement'', or ''complement_object''.', opt.IRMode);
end

% ---------------- Expand 2D masks to 3D ----------------
maskROI = repmat(roi2D, [1 1 Nz]);
maskIR  = repmat(ir2D,  [1 1 Nz]);

% Remove overlap if any
if any(maskROI(:) & maskIR(:))
    warning('computeROVirFromEyeROI:MaskOverlap', 'ROI and IR overlap. Removing overlap from IR.');
    maskIR(maskROI) = false;
end

idxO = maskROI(:);
idxG = maskIR(:);

if nnz(idxO) == 0, error('computeROVirFromEyeROI:EmptyROI', 'ROI mask has 0 voxels.'); end
if nnz(idxG) == 0, error('computeROVirFromEyeROI:EmptyIR',  'IR mask has 0 voxels.'); end

% ---------------- Form A and B ----------------
G = reshape(g_cal, [], nCh);   % [Nvox, nCh]
G_O = G(idxO, :);
G_G = G(idxG, :);

A = (G_O' * G_O) / nnz(idxO);
B = (G_G' * G_G) / nnz(idxG);

% enforce Hermitian
A = (A + A') / 2;
B = (B + B') / 2;
if debugFlag
    % Diagnostics
    fprintf('Nroi=%d, Nir=%d\n', nnz(idxO), nnz(idxG));
    fprintf('trace(A)=%.3e, trace(B)=%.3e\n', real(trace(A)), real(trace(B)));
    fprintf('rcond(B)=%.3e\n', rcond(B));
    
    eA = eig(A); eB = eig(B);
    fprintf('eig(A) min/median/max = %.3e / %.3e / %.3e\n', min(real(eA)), median(real(eA)), max(real(eA)));
    fprintf('eig(B) min/median/max = %.3e / %.3e / %.3e\n', min(real(eB)), median(real(eB)), max(real(eB)));
end
% ---------------- Regularize B ----------------
tB = real(trace(B));
if ~isfinite(tB) || tB <= 0
    tB = norm(B, 'fro');
end
epsReg = opt.EpsFactor * (tB / nCh);
if ~isfinite(epsReg) || epsReg == 0
    epsReg = opt.EpsFactor;
end
Breg = B + epsReg * eye(nCh);

% ---------------- GEVD ----------------
[W, d] = eig(A, Breg, 'vector');
lambda = real(d);
[lambda, order] = sort(lambda, 'descend');
W = W(:, order);
W = W ./ vecnorm(W, 2, 1);

% ---------------- Return masks (optional) ----------------
masksOut = struct();
masksOut.roi2D = roi2D;
masksOut.ir2D  = ir2D;
masksOut.roi3D = maskROI;
masksOut.ir3D  = maskIR;

end

