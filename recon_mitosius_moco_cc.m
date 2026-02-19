function [mDir, C_out, ccInfo] = recon_mitosius_moco_cc( ...
    RECON_FOLDER_DATE, meas_name, seqName, CfilePath, woBin, maskFilePath, varargin)
% recon_mitosius_moco_cc
%
% Purpose (HPC):
%   - Load raw data + pulseq trajectory
%   - Normalize y_tot
%   - Apply rigid motion correction (rotation on t, translation via phase ramp on y)
%   - Prepare eye mask
%   - OPTIONAL coil compression for faster bmMitosius_create downstream
%   - Export mitosius folder
%
% I/O (similar to recon_mitosius):
%   Inputs:
%     RECON_FOLDER_DATE : e.g. '260201'
%     meas_name         : raw data file name, contains MIDxxxx
%     seqName           : pulseq .seq file name
%     CfilePath         : path to C.mat (must contain variable C)
%     woBin             : 1/0 (default 1)
%     maskFilePath      : if woBin==0, path to mask file (no extension needed is OK)
%
%   Outputs:
%     mDir   : output mitosius directory
%     C_out  : (compressed or original) sensitivity maps
%     ccInfo : coil compression diagnostics (empty if compression off)
%
% Name-Value Options (motion correction):
%   'doMotionCorrection' (false)
%   'rpFilePath'         ('')   : path to rp_*.txt (6 columns: Tx Ty Tz Rx Ry Rz)
%   'seqBinningMatPath'  ('')   : path to sequentialBinning_winXXs.mat (contains variable 'mask')
%   'doRot'              (true)
%   'doTrans'            (true)
%   'rotSigns'           ([1 1 -1])
%   'tranSigns'          ([-1 1 1])
%   'nShotOff'           (14)   : overrides reader.acquisitionParams.nShot_off
%
% Name-Value Options (coil compression):
%   'doCoilCompression' (false)
%   'nChCompressed'     (20)
%   'ccMethod'          (2)    : 0=SVD, 1=energy, 2=ROI, 3=topN mat
%   'weightsNormPath'    ('')   : required if ccMethod=2
%   'saveDiagnostics'   (true)
%
% Notes:
%   - This refactors your S02 rigid motion correction logic.
%   - Motion correction uses the sequential-binning mask to assign one rigid transform per line.
%   - Trajectory rotation is applied line-wise; translation is applied as a phase ramp on y.
%   - Coil compression is applied *after* motion correction, before bmMitosis/bmMitosius_create.

%% ---------------------- Defaults ----------------------
if nargin < 5 || isempty(woBin)
    woBin = 1;
end
if nargin < 6
    maskFilePath = '';
end

% Parse NV options with inputParser (HPC-safe)
ip = inputParser;
ip.KeepUnmatched = false;

% Motion correction options
addParameter(ip, 'doMotionCorrection', false, @(x)islogical(x) || isnumeric(x));
addParameter(ip, 'rpFilePath', '', @(x)ischar(x) || isstring(x));
addParameter(ip, 'seqBinningMatPath', '', @(x)ischar(x) || isstring(x));
addParameter(ip, 'doRot', true, @(x)islogical(x) || isnumeric(x));
addParameter(ip, 'doTrans', true, @(x)islogical(x) || isnumeric(x));
addParameter(ip, 'rotSigns', [1 1 -1], @(x)isnumeric(x) && numel(x)==3);
addParameter(ip, 'tranSigns', [-1 1 1], @(x)isnumeric(x) && numel(x)==3);
addParameter(ip, 'nShotOff', 14, @(x)isnumeric(x) && isscalar(x));

% Coil compression options
addParameter(ip, 'doCoilCompression', false, @(x)islogical(x) || isnumeric(x));
addParameter(ip, 'nChCompressed', 20, @(x)isnumeric(x) && isscalar(x) && x>0);
addParameter(ip, 'ccMethod', 2, @(x)isnumeric(x) && isscalar(x));
addParameter(ip, 'weightsNormPath', '', @(x)ischar(x) || isstring(x));
addParameter(ip, 'saveDiagnostics', true, @(x)islogical(x) || isnumeric(x));

parse(ip, varargin{:});
opt = ip.Results;

%% ---------------------- Print mode ----------------------
disp('=================================================')
disp('recon_mitosius_moco (HPC)')
disp(['meas_name: ', char(meas_name)])
disp(['seqName:   ', char(seqName)])

if woBin
    disp('Mask mode: woBin (auto-generated mask)')
else
    disp('Mask mode: using input maskFilePath:')
    disp(char(maskFilePath))
end

%% ---------------------- Initialize directories ----------------------
token = regexp(meas_name, 'MID\d+', 'match');
measMID = token{1};
disp(['measMID: ', measMID])

datasetDir = ['/usr/src/app/datasets/' RECON_FOLDER_DATE '/'];
seqFolder  = ['/usr/src/app/datasets/' RECON_FOLDER_DATE '/'];
reconDir   = ['/usr/src/app/recon_folder/' RECON_FOLDER_DATE '/' measMID '_recon'];

disp(['reconDir: ', reconDir])

if ~isfolder(reconDir)
    mkdir(reconDir);
    disp(['Directory created: ', reconDir]);
end
cd(reconDir);

% otherDir layout
if woBin
    otherDir = [reconDir, '/T1_LIBRE_woBinning/other/'];
else
    otherDir = [reconDir, '/moco/other/'];
end
if ~isfolder(otherDir)
    mkdir(otherDir);
    disp(['Directory created: ', otherDir]);
end

measureFile = fullfile(datasetDir, meas_name);
seqFile     = fullfile(seqFolder,  seqName);

disp(['measureFile: ', measureFile])
disp(['seqFile:     ', seqFile])

%% ---------------------- Load and configure data ----------------------
seqParams = extract_seq_params(seqFile);

reader = createRawDataReader(measureFile, 1);

% Acquisition from Bern
reader.acquisitionParams.nShot_off = opt.nShotOff;
reader.acquisitionParams.traj_type = 'pulseq';
reader.acquisitionParams.pulseqTrajFile_name = char(seqFile);

% Optional consistency check
try
    isMatch = check_hash(measureFile, reader.acquisitionParams.pulseqTrajFile_name); %#ok<NASGU>
catch
    disp('check_hash failed or not found; continuing.')
end

if isfield(seqParams, 'nshot'); reader.acquisitionParams.nShot = seqParams.nshot; end
if isfield(seqParams, 'nseg');  reader.acquisitionParams.nSeg  = seqParams.nseg;  end

% Read raw (filters nShotOff + SI as in your code)
y_tot = reader.readRawData(true, true);

% Trajectory
t_tot  = bmTraj(reader.acquisitionParams);
ve_tot = bmVolumeElement(t_tot, 'voronoi_full_radial3'); %#ok<NASGU>

% Load sensitivity maps
S = load(CfilePath, 'C');
C = S.C;
disp(['C loaded from: ', CfilePath])

%% ---------------------- Normalize (same as your HPC recon_mitosius) ----------------------
matrix_size = 120;
N_u = [matrix_size, matrix_size, matrix_size];

% Your dK_u formula used seqParams.fov*2e3
dK_u = 1 ./ (seqParams.fov * 2e3);

C_preview = bmImResize(C, [48, 48, 48], N_u);

disp('Computing bmMathilda preview for normalization...')
x_tot = bmMathilda(y_tot, t_tot, ve_tot, C_preview, N_u, N_u, dK_u);

temp_im = x_tot( ...
    round(matrix_size/4):round(matrix_size/4*3), ...
    round(matrix_size/4):round(matrix_size/4*3), ...
    round(matrix_size/2));

normalize_val = mean(temp_im(:));
disp(['normalize_val: ', num2str(normalize_val)])

if mean(abs(y_tot(:))) < 1
    y_tot = y_tot / normalize_val;
    disp('y_tot normalized by normalize_val.')
else
    disp('y_tot seems already normalized (real(y_tot) >= 1).')
end

%% ---------------------- Motion correction (rigid) ----------------------
t_tot_clean = t_tot;
y_tot_clean = y_tot;

if logical(opt.doMotionCorrection)
    disp('-------------------------------------------------')
    disp('Rigid motion correction: ENABLED')

    % If user did not pass paths, we try common HPC layout guesses
    rpFilePath = char(opt.rpFilePath);
    seqBinningMatPath = char(opt.seqBinningMatPath);

    if isempty(rpFilePath)
        % Try: reconDir/moco/output/rp_*.txt
        tryHit = dir(fullfile(reconDir, 'moco', 'output', 'rp_*.txt'));
        if ~isempty(tryHit)
            rpFilePath = fullfile(tryHit(1).folder, tryHit(1).name);
            disp(['Auto-found rpFilePath: ', rpFilePath])
        end
    end

    if isempty(seqBinningMatPath)
        % Try: reconDir/moco/sequentialBinning_*.mat
        tryHit = dir(fullfile(reconDir, 'moco', 'sequentialBinning_*.mat'));
        if ~isempty(tryHit)
            seqBinningMatPath = fullfile(tryHit(1).folder, tryHit(1).name);
            disp(['Auto-found seqBinningMatPath: ', seqBinningMatPath])
        end
    end

    [t_tot_clean, y_tot_clean] = apply_rigid_moco( ...
        reader, t_tot, y_tot, seqParams, rpFilePath, seqBinningMatPath, ...
        logical(opt.doRot), logical(opt.doTrans), opt.rotSigns, opt.tranSigns);

else
    disp('-------------------------------------------------')
    disp('Rigid motion correction: DISABLED')
end

%% ---------------------- Mask path (woBin or provided) ----------------------
if woBin == 1
    disp('Generating woBin eMask...')
    eMask = ones(1, seqParams.nshot * seqParams.nseg) > 0;
    eMask(1:reader.acquisitionParams.nShot_off * seqParams.nseg) = 0;

    eMaskSavePath = fullfile(otherDir, 'eMask_woBin.mat');
    save(eMaskSavePath, 'eMask');
    disp(['eMask saved: ', eMaskSavePath])

    % For loading, allow without extension (your style)
    allLinesBinningPath = fullfile(otherDir, 'eMask_woBin');
else
    disp('Using provided maskFilePath:')
    disp(char(maskFilePath))
    allLinesBinningPath = char(maskFilePath);
end

%% ---------------------- Mitosius output folder ----------------------
if woBin
    mDir = [reconDir, '/T1_LIBRE_woBinning/mitosius/'];
else
    mDir = [reconDir, '/moco/mitosius/'];
end

if ~isfolder(mDir)
    mkdir(mDir);
else
    mDir0 = char(string(mDir));           % force to scalar text, then char row
    mDir0 = strtrim(mDir0);               % remove whitespace
    mDir = regexprep(mDir0, '[\\/]+$', '');   % remove trailing / or \
    mDirParent = fileparts(mDir);
    mDir = char(fullfile(mDirParent, ['mitosius_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]));
    mkdir(mDir);
end
disp(['mDir: ', mDir])

%% ---------------------- Write Mitosius (with optional coil compression) ----------------------
[C_out, ccInfo] = write_mitosius_from_mask( ...
    t_tot_clean, y_tot_clean, allLinesBinningPath, reader.acquisitionParams, C, mDir, ...
    'doCoilCompression', logical(opt.doCoilCompression), ...
    'nChCompressed',     opt.nChCompressed, ...
    'ccMethod',          opt.ccMethod, ...
    'weightsNormPath',   char(opt.weightsNormPath), ...
    'saveDiagnostics',   logical(opt.saveDiagnostics) );

disp('=================================================')
disp('Done.')
disp(['Mitosius saved at: ', mDir])

end


%% ========================================================================
% Helper: apply rigid motion correction (refactor of S02_motion_correction_volume_wise)
function [t_tot_clean, y_tot_clean] = apply_rigid_moco( ...
    reader, t_tot, y_tot, seqParams, rpFilePath, seqBinningMatPath, doRot, doTrans, rotSigns, tranSigns)

% Validate input paths
if isempty(rpFilePath) || ~exist(rpFilePath, 'file')
    error('rpFilePath missing or not found: %s', rpFilePath);
end
if isempty(seqBinningMatPath) || ~exist(seqBinningMatPath, 'file')
    error('seqBinningMatPath missing or not found: %s', seqBinningMatPath);
end

disp(['rpFilePath: ', rpFilePath])
disp(['seqBinningMatPath: ', seqBinningMatPath])

% Read rp table (Tx Ty Tz Rx Ry Rz)
rpTable = importdata(rpFilePath);
TransParam = [rpTable(:,1)'; rpTable(:,2)'; rpTable(:,3)'];
RotParam   = [rpTable(:,4)'; rpTable(:,5)'; rpTable(:,6)'];

% Load sequential binning mask (expects variable "mask" in .mat)
M = load(seqBinningMatPath);
if ~isfield(M, 'mask')
    error('sequential binning mat does not contain variable "mask": %s', seqBinningMatPath);
end
mask = M.mask;

acq = reader.acquisitionParams;

nSeg     = acq.nSeg;
nShot    = acq.nShot;
nShotOff = acq.nShot_off;
nLines   = acq.nLine;
nCol     = acq.N;

% --- Expand motion params to per-line using mask weights (as in your S02 code)
RotParam_volume   = repmat(RotParam,   [1, 1, nLines]);
TransParam_volume = repmat(TransParam, [1, 1, nLines]);

% mask: [nbins, nSeg*nShot] or similar -> S02 assumes reshapeable to [nbins,nSeg,nShot]
size_Mask = size(mask);
nbins = size_Mask(1);

mask_temp   = reshape(mask, [nbins, nSeg, nShot]);
mask_temp   = reshape(mask_temp, [1, size(mask_temp,1), size(mask_temp,2)*size(mask_temp,3)]);
mask_expand = repmat(mask_temp, [3, 1, 1]);

% Weighted sum across bins
RotParam_volume = RotParam_volume .* mask_expand;
RotParam_volume = sum(RotParam_volume, 2);
RotParam_volume = squeeze(RotParam_volume);

TransParam_volume = TransParam_volume .* mask_expand;
TransParam_volume = sum(TransParam_volume, 2);
TransParam_volume = squeeze(TransParam_volume);

% Reshape to [3, nSeg, nShot] then remove shotOff and SI segment
RotParam_volume = reshape(RotParam_volume, 3, nSeg, nShot);
RotParam_volume(:, :, 1:nShotOff) = [];
RotParam_volume(:, 1, :) = [];        % remove SI segment
RotParam_volume = reshape(RotParam_volume, 3, []);

TransParam_volume = reshape(TransParam_volume, 3, nSeg, nShot);
TransParam_volume(:, :, 1:nShotOff) = [];
TransParam_volume(:, 1, :) = [];      % remove SI segment
TransParam_volume = reshape(TransParam_volume, 3, []);

% --- Rotation correction: rotate trajectory line-wise
if doRot
    disp('Motion correction: ROTATION ON')
    t_tot_clean = zeros(size(t_tot), 'like', t_tot);

    for iLine = 1:size(RotParam_volume,2)
        Rot_mat = euler_to_rotation(RotParam_volume(:, iLine), rotSigns);

        te = reshape(t_tot(:,:,iLine), 3, []);
        te_rot = Rot_mat * te;
        te_rot = reshape(te_rot, 3, nCol, []);

        t_tot_clean(:,:,iLine) = te_rot;
    end
else
    disp('Motion correction: ROTATION OFF')
    t_tot_clean = t_tot;
end

% --- Translation correction: phase ramp on y_tot (uses ORIGINAL t_tot like your S02)
if doTrans
    disp('Motion correction: TRANSLATION ON')

    kx = t_tot(1,:,:);
    ky = t_tot(2,:,:);
    kz = t_tot(3,:,:);

    tempx = kx .* reshape(tranSigns(1) * TransParam_volume(1,:), 1,1,[]);
    tempy = ky .* reshape(tranSigns(2) * TransParam_volume(2,:), 1,1,[]);
    tempz = kz .* reshape(tranSigns(3) * TransParam_volume(3,:), 1,1,[]);

    PhaseOffset = exp(-2*pi*1i*(tempx + tempy + tempz));
    y_tot_clean = y_tot .* PhaseOffset;

else
    disp('Motion correction: TRANSLATION OFF')
    y_tot_clean = y_tot;
end

disp('Rigid motion correction applied.')

end


%% ========================================================================
% Helper: write mitosius from mask, with optional coil compression
function [C_out, ccInfo] = write_mitosius_from_mask( ...
    t_tot_clean, y_tot_clean, allLinesBinningPath, acquisitionParams, C, mDir, varargin)

% Parse NV options
ip = inputParser;
ip.KeepUnmatched = false;

addParameter(ip, 'doCoilCompression', false, @(x)islogical(x) || isnumeric(x));
addParameter(ip, 'nChCompressed', 20, @(x)isnumeric(x) && isscalar(x) && x>0);
addParameter(ip, 'ccMethod', 0, @(x)isnumeric(x) && isscalar(x));
addParameter(ip, 'weightsNormPath', '', @(x)ischar(x) || isstring(x));
addParameter(ip, 'saveDiagnostics', true, @(x)islogical(x) || isnumeric(x));

parse(ip, varargin{:});
opt = ip.Results;

t_tot = t_tot_clean;
y_tot = y_tot_clean;
p = acquisitionParams;

% Load mask
tempMask = load(allLinesBinningPath);
fields = fieldnames(tempMask);
tempMask = tempMask.(fields{1});

disp(['Mask loaded: ', char(allLinesBinningPath)])

% Prepare mask to match filtered data (remove SI + shotOff)
nbins = size(tempMask, 1);
tempMask = reshape(tempMask, [nbins, p.nSeg, p.nShot]);
tempMask(:, 1, :) = [];                 % remove SI segment
tempMask(:, :, 1:p.nShot_off) = [];     % remove shots-off
tempMask = bmPointReshape(tempMask);

% Optional coil compression
C_out = C;
ccInfo = struct();

if logical(opt.doCoilCompression)
    disp('Coil compression: ENABLED')

    [y_tot, C_out, ccInfo] = coilCompressCore( ...
        y_tot, C_out, ...
        'method',          opt.ccMethod, ...
        'nChCompressed',   opt.nChCompressed, ...
        'weightsNormPath', char(opt.weightsNormPath), ...
        'saveDiagnostics', logical(opt.saveDiagnostics), ...
        'mDir',            mDir);

    disp(['Compressed coils: y_tot -> ', num2str(size(y_tot,1)), ...
          ', C_out -> ', num2str(size(C_out,4))])
    
    CoutfilePath = char(fullfile(mDir, sprintf('C_compressed_%d.mat', opt.nChCompressed)));
    save(CoutfilePath, 'C_out', '-v7.3'); % save C_out 
    disp(['C_out has been saved here: ', CoutfilePath]);
else
    disp('Coil compression: DISABLED')
end

% Run mitosis and export mitosius
[y, t] = bmMitosis(y_tot, t_tot, tempMask);
y = bmPermuteToCol(y);
ve = bmVolumeElement(t, 'voronoi_full_radial3');

bmMitosius_create(mDir, y, t, ve);

disp('Mitosius files saved:')
disp(mDir)

end




%% ========================================================================
% Helper: Euler -> rotation matrix (same as your S02)
function R = euler_to_rotation(R_array, rotSigns)
Rx = rotSigns(1) * R_array(1);
Ry = rotSigns(2) * R_array(2);
Rz = rotSigns(3) * R_array(3);

Rx_mat = [1 0 0;
          0 cos(Rx) -sin(Rx);
          0 sin(Rx)  cos(Rx)];

Ry_mat = [cos(Ry) 0 sin(Ry);
          0 1 0;
         -sin(Ry) 0 cos(Ry)];

Rz_mat = [cos(Rz) -sin(Rz) 0;
          sin(Rz)  cos(Rz) 0;
          0 0 1];

R = Rz_mat * Ry_mat * Rx_mat;
R = R.'; % inverse to correct back
end