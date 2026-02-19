%% ========================================================================
% Helper: coil compression core (SVD / energy / ROI-selection / GEVD)
function [y_tot_comp, C_comp, info] = coilCompressCore(y_tot, C, varargin)

ip = inputParser;
ip.KeepUnmatched = false;

addParameter(ip, 'method', 0, @(x)isnumeric(x) && isscalar(x));
addParameter(ip, 'nChCompressed', 20, @(x)isnumeric(x) && isscalar(x) && x>0);
addParameter(ip, 'weightsNormPath', '', @(x)ischar(x) || isstring(x));
addParameter(ip, 'saveDiagnostics', true, @(x)islogical(x) || isnumeric(x));
addParameter(ip, 'mDir', '', @(x)ischar(x) || isstring(x));
parse(ip, varargin{:});
opt = ip.Results;

info = struct();
info.method = opt.method;
info.nChCompressed = opt.nChCompressed;

nCh = size(y_tot, 1);
if size(C, 4) ~= nCh
    error('Mismatch: y_tot has %d coils, but C has %d coils.', nCh, size(C,4));
end

switch opt.method
    case 0
        % ---------- SVD virtual coils ----------
        % nx = size(y_tot, 2);
        disp('Coil compression method: SVD virtual coils');
        % ntviews = size(y_tot, 3);
        nCh     = size(y_tot, 1);
        nSample = size(y_tot, 2);
        nLine   = size(y_tot, 3);

        fprintf('nSample: %d, nCh: %d, nLine: %d\n', nSample, nCh, nLine);
        D = reshape(permute(y_tot, [2 3 1]), [], nCh);   % (nSample*nLine) x nCh, each column = one coil

        [~, S, V] = svd(D, 'econ');
        k = opt.nChCompressed;
        V_k = V(:, 1:k);
        % Compress data
        Dcomp = D * V_k;                                  % [(nSample*nLine), k]
        y_tot_comp = permute(reshape(Dcomp, nSample, nLine, k), [3 1 2]);  % -> [k, nSample, nLine]

        % Compress C
        Cmat = reshape(C, [], nCh);            % [Nvox, nCh]
        Cmat_comp = Cmat * V_k;                % [Nvox, k]
        C_comp = reshape(Cmat_comp, size(C,1), size(C,2), size(C,3), k);

        % Diagnostics: explained variance
        singular_values = diag(S);
        total = sum(singular_values.^2);
        percentage_expl = (singular_values(1:k).^2 / total) * 100;
        info.percentage_explained = percentage_expl;
        info.V = V;


        if logical(opt.saveDiagnostics) && ~isempty(char(opt.mDir))
            f = figure;
            plot(percentage_expl, '.-', 'LineWidth', 2, 'MarkerSize', 14);
            xlabel('Virtual coil #');
            ylabel('Explained variance [%]');
            title(sprintf('SVD coil compression (N=%d), sum=%.2f%%', opt.nChCompressed, sum(percentage_expl)));
            outPng = fullfile(char(opt.mDir), sprintf('coilCompression_svd_expVar_%d.png', opt.nChCompressed));
            saveas(f, outPng);
            close(f);
            info.expVarFig = outPng;
        end

    case 1
        % ---------- Energy-based top-N selection ----------
        coilEnergy = squeeze(sum(abs(y_tot).^2, [2 3]));
        [~, idxSorted] = sort(coilEnergy, 'descend');
        idxKeep = idxSorted(1:opt.nChCompressed);

        y_tot_comp = y_tot(idxKeep, :, :);
        C_comp = C(:,:,:, idxKeep);

        info.idxKeep = idxKeep;
        info.coilEnergy = coilEnergy;

    case 2
        % ---------- ROI-weight top-N selection ----------
        weightsNormPath = char(opt.weightsNormPath);
        if isempty(weightsNormPath)
            error('ccMethod=2 requires weightsNormPath for coil ranking.');
        end

        % robust call: try (reconDir,C), fallback to (reconDir)
        w_cell = load(weightsNormPath);
        weights_norm = w_cell.w;

        [~, idxSorted] = sort(weights_norm, 'descend');
        idxKeep = idxSorted(1:opt.nChCompressed);

        y_tot_comp = y_tot(idxKeep, :, :);
        C_comp = C(:,:,:, idxKeep);

        info.idxKeep = idxKeep;
        info.weights_norm = weights_norm;

        if logical(opt.saveDiagnostics) && ~isempty(char(opt.mDir))
            baseOut = fileparts(char(opt.mDir));
            save(fullfile(baseOut, 'weights_norm.mat'), 'weights_norm', '-v7.3');
            save(fullfile(baseOut, 'idx_coilSelection.mat'), 'idxSorted', '-v7.3');
        end

    case 3
     % ---------- ROVir / GEVD virtual coils ----------
        disp('Coil compression method: ROVir (GEVD) virtual coils');
        weightsNormPath = char(opt.weightsNormPath);
        if isempty(weightsNormPath)
            error('ccMethod=3 requires weightsNormPath for ROVir.');
        end
        nCh     = size(y_tot, 1);
        nSample = size(y_tot, 2);
        nLine   = size(y_tot, 3);

        % robust call: try (reconDir,C), fallback to (reconDir)
        w_cell = load(weightsNormPath);
        W = w_cell.W;
        % Pick Nv
        Nv = opt.nChCompressed;                    % simplest: keep top-Nv
        Wk = W(:, 1:Nv);                           % [nCh, Nv]

        % ---- Apply to raw data (same reshaping pattern as your SVD code) ----
        Draw  = reshape(permute(y_tot, [2 3 1]), [], nCh);   % [(nSample*nLine), nCh]
        [Wk_ortho,~] = qr(Wk,0);
        Dcomp     = Draw * conj(Wk_ortho);

        y_tot_comp = permute(reshape(Dcomp, nSample, nLine, Nv), [3 1 2]);
        
        % ---- Apply same mixing to coil sensitivities / C  ----
        Cmat = reshape(C, [], nCh);
        % Cmat_comp = Cmat * Wk;
        Cmat_comp = Cmat * conj(Wk_ortho);
        C_comp = reshape(Cmat_comp, size(C,1), size(C,2), size(C,3), Nv);
        
        % ---- Diagnostics: per-virtual-coil SIR values ----
        info.weights_norm = W;


    otherwise
        error('Unknown ccMethod: %d', opt.method);
end

end
