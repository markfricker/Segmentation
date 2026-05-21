%%
% sweepCellpose.m
%
% Sweep Cellpose hyperparameters to identify optimal settings for
% mitochondria / organelle segmentation.
%
% Unlike sweepSegmentation (which treats Cellpose as fixed GT and sweeps
% classical algorithms), this script is entirely about tuning Cellpose
% itself.  There is no external ground truth — the output metrics
% (object count, foreground coverage) and visual inspection are used
% to identify the 'Goldilocks' parameter set: not under-segmenting
% faint/small objects, but not over-merging or producing noise detections.
%
% IMAGE INPUT
%   Priority order:
%     1. Workspace variable named by 'workspaceVar' (set to '' to skip).
%        Any 2-D or 2-D single-channel numeric array is accepted.
%     2. File at blobImagePath/blobImageFile.
%   If neither is available the script errors with a clear message.
%
% CELLPOSE STRATEGY
%   For each (model, nIter) combination:
%     Sweep CellThreshold (cellProb) × FlowErrorThreshold (flowThreshold).
%     Record n_objects and foreground coverage for each combo.
%   Output a 2-D console table per model, plus heatmap and visual-grid
%   figures for fast identification of the best parameter region.
%
% BLOB FALLBACK
%   If the Cellpose add-on or cellposeSegment.m is not available, the
%   script falls through to a classical blob segmentation sweep using
%   watershedSegment (marker-controlled watershed on the enhanced image).
%   Sweeps threshold × hMinima in the same display format as the Cellpose
%   sweep, so results are directly comparable when Cellpose becomes available.
%
% CACHING
%   Cellpose sweep results (Lall, nObjAll, covAll) are cached in the
%   workspace.  Re-running the script with the same sweep grid but
%   different display options skips the Cellpose calls entirely.
%   Blob results (blobLall, blobNObjAll, blobCovAll) are cached separately.
%   To force a full rerun, clear the workspace or change a grid parameter.
%
% METRICS
%   n objects   — number of distinct detected objects (= max label index).
%   coverage %  — fraction of image pixels labelled as foreground (× 100).
%
% FIGURES (per model × nIter combination)
%   Fig A — Combined heatmaps: n_objects (left) and coverage % (right).
%   Fig B — Visual label grid (full image).
%
% CONFIGURATION
%   cpModels  — cell array of Cellpose model names to test.
%   cpDiam    — object diameter in pixels (fixed; tune separately).
%   cpVals    — CellThreshold (cellProb) sweep values.  Range: [-6, 6].
%   ftVals    — FlowErrorThreshold sweep values.  Range: [0.1, 3].
%   niVals    — nIter values.  0 = Cellpose default (~200).
%
% REQUIREMENTS
%   cellposeSegment.m on path (Segmentation_sandbox/src).
%   watershedSegment.m on path (Segmentation_sandbox/src) for blob fallback.
%   Medical Imaging Toolbox Interface for Cellpose Library Add-On.
%   Python environment with cellpose installed.

clc; close all;   % intentionally no 'clear' — workspace cache is preserved

% =========================================================================
% CONFIGURATION
% =========================================================================
blobFunctionPath = 'C:\Users\dops0035\Documents\Research\Matlab Projects\BlobFilters_sandbox\src';
blobImagePath    = 'C:\Users\dops0035\Documents\Research\Matlab Projects\BlobFilters_sandbox\demos';
blobImageFile    = 'mitImagecrop.mat';

% ---- Workspace image override -------------------------------------------
% Set to a variable name (string) that already exists in the base workspace
% (e.g. created by running the AnalyzER organelle enhance step).
% Set to '' to always load from blobImagePath/blobImageFile.
workspaceVar = 'im2D';

% ---- Cellpose models and nIter to test ----------------------------------
cpModels  = {'cyto3', 'bact_fluor_cp3'};
niVals    = [0, 2000];
cpDiam    = 10;           % expected object diameter in pixels (fixed)

% Combinations to exclude: each entry is {modelName, nIterValue}.
skipCombos = {{'bact_fluor_cp3', 2000}};

% ---- Cellpose sweep grid ------------------------------------------------
cpVals    = [-1, 0, 1];          % CellThreshold (cellProb) sweep values
ftVals    = [0.60, 0.80, 1.00]; % FlowErrorThreshold sweep values
cpMinSize = 64;                  % px^2 — discard objects below this area

% ---- Blob fallback sweep grid -------------------------------------------
% Used when Cellpose is not available.  Rows = threshold, Cols = hMinima.
blobThreshVals = [0.15, 0.25, 0.35];  % foreground threshold (image [0,1])
blobHMinVals   = [1, 3, 5];           % h-minima suppression depth
blobMinSize    = 20;                   % px^2 min object area

% ---- Figure export -------------------------------------------------------
doExport = false;
outDir   = fullfile(blobFunctionPath, '..', 'docs', 'figures');

% =========================================================================
% Path setup
% =========================================================================
demoDir = fileparts(mfilename('fullpath'));
addpath(fullfile(demoDir, '..', 'src'));
if ~isempty(blobFunctionPath) && exist(blobFunctionPath, 'dir')
    addpath(blobFunctionPath);
end

% Detect what's available
hasCellpose = exist('cellposeSegment', 'file') == 2 && ...
              exist('cellpose',        'file') ~= 0;
hasWatershed = exist('watershedSegment', 'file') == 2;

if ~hasCellpose
    fprintf('=== Cellpose not available — running blob fallback sweep ===\n\n');
    if ~hasWatershed
        error('sweepCellpose:noFallback', ...
              'Neither cellposeSegment nor watershedSegment found on path.\n%s\n%s', ...
              'Cellpose: install "Medical Imaging Toolbox Interface for Cellpose Library"', ...
              'Blob fallback: add Segmentation_sandbox/src to path.');
    end
end

% =========================================================================
% 1.  Load image
% =========================================================================
Ireal = loadImage(workspaceVar, blobImagePath, blobImageFile);

% =========================================================================
% Branch: Cellpose sweep  OR  blob fallback
% =========================================================================
if hasCellpose
    runCelloseSweep();
else
    runBlobFallback();
end

fprintf('\nDone.\n');

% =========================================================================
% NESTED FUNCTIONS  (share variables with the script workspace)
% =========================================================================

% -------------------------------------------------------------------------
function runCelloseSweep()
% Full Cellpose parameter sweep — mirrors the original sweepCellpose logic.

    nM  = numel(cpModels);
    nNI = numel(niVals);
    nCP = numel(cpVals);
    nFT = numel(ftVals);

    % Build workspace cache tag
    cpCacheTagNew = sprintf('%s|ni%s|d%.0f|cp%s|ft%s|ms%.0f', ...
        strjoin(cpModels, '_'), num2str(niVals,'%g '), cpDiam, ...
        num2str(cpVals,'%.2f '), num2str(ftVals,'%.2f '), cpMinSize);

    doSweep = true;
    if evalin('base','exist(''cpCacheTag'',''var'')') && ...
       strcmp(evalin('base','cpCacheTag'), cpCacheTagNew) && ...
       evalin('base','exist(''Lall'',''var'')') && ...
       isequal(size(evalin('base','Lall')), [nM nNI nCP nFT])
        fprintf('=== Using cached Cellpose results ===\n');
        fprintf('    Change a grid parameter or clear Lall to force rerun.\n\n');
        doSweep = false;
        Lall    = evalin('base','Lall');
        nObjAll = evalin('base','nObjAll');
        covAll  = evalin('base','covAll');
    end

    if doSweep
        nObjAll = nan(nM, nNI, nCP, nFT);
        covAll  = nan(nM, nNI, nCP, nFT);
        Lall    = cell(nM, nNI, nCP, nFT);

        nTotal = nM * nNI * nCP * nFT;
        fprintf('\n=== Cellpose parameter sweep  (%d configurations) ===\n', nTotal);

        for mi = 1:nM
            for ni = 1:nNI
                niStr = '';
                if niVals(ni) > 0, niStr = sprintf('  nIter=%d', niVals(ni)); end

                if any(cellfun(@(s) strcmp(cpModels{mi},s{1}) && niVals(ni)==s{2}, skipCombos))
                    fprintf('\n--- %s%s --- [SKIPPED]\n', cpModels{mi}, niStr);
                    continue;
                end

                fprintf('\n--- %s%s ---\n', cpModels{mi}, niStr);
                t0 = tic;

                for ci = 1:nCP
                    for fi = 1:nFT
                        pCP.model         = cpModels{mi};
                        pCP.diameter      = cpDiam;
                        pCP.cellProb      = cpVals(ci);
                        pCP.flowThreshold = ftVals(fi);
                        pCP.nIter         = niVals(ni);
                        pCP.minSize       = cpMinSize;

                        fprintf('  cp=%+d  ft=%.2f  ... ', cpVals(ci), ftVals(fi));
                        tic;
                        try
                            [Lk, ~] = cellposeSegment(Ireal, pCP);
                            n_k = double(max(Lk(:)));
                            c_k = 100 * nnz(Lk > 0) / numel(Lk);
                            nObjAll(mi,ni,ci,fi) = n_k;
                            covAll(mi,ni,ci,fi)  = c_k;
                            Lall{mi,ni,ci,fi}    = Lk;
                            fprintf('%.1fs  n=%3d  cov=%5.1f%%\n', toc, n_k, c_k);
                        catch ME
                            fprintf('FAILED (%s)\n', ME.message);
                            Lall{mi,ni,ci,fi} = zeros(size(Ireal),'uint16');
                        end
                    end
                end

                fprintf('  Done (%.1fs total)\n', toc(t0));
            end
        end

        % Store in base workspace for caching
        assignin('base', 'Lall',       Lall);
        assignin('base', 'nObjAll',    nObjAll);
        assignin('base', 'covAll',     covAll);
        assignin('base', 'cpCacheTag', cpCacheTagNew);
    end

    % ---- Console tables -------------------------------------------------
    for mi = 1:nM
        for ni = 1:nNI
            if any(cellfun(@(s) strcmp(cpModels{mi},s{1}) && niVals(ni)==s{2}, skipCombos))
                continue;
            end
            niStr = '';
            if niVals(ni) > 0, niStr = sprintf('  nIter=%d', niVals(ni)); end
            printCPtable(cpModels{mi}, niStr, cpVals, ftVals, ...
                         squeeze(nObjAll(mi,ni,:,:)), ...
                         squeeze(covAll(mi,ni,:,:)));
        end
    end

    % ---- Figures --------------------------------------------------------
    [imgH, imgW] = size(Ireal);
    figNum = 0;

    for mi = 1:nM
        for ni = 1:nNI
            if any(cellfun(@(s) strcmp(cpModels{mi},s{1}) && niVals(ni)==s{2}, skipCombos))
                continue;
            end
            niStr = '';
            if niVals(ni) > 0, niStr = sprintf(' nIter=%d', niVals(ni)); end
            mTag = sprintf('%s%s', cpModels{mi}, niStr);

            nObj_mn = squeeze(nObjAll(mi,ni,:,:));
            cov_mn  = squeeze(covAll(mi,ni,:,:));
            L_mn    = squeeze(Lall(mi,ni,:,:));

            % Heatmap figure
            figNum = figNum + 1;
            hfa = figure(figNum);
            set(hfa, 'Name', sprintf('Cellpose heatmaps — %s', mTag), ...
                     'NumberTitle','off', 'Color',[0.13 0.13 0.13]);
            tla = tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','compact');
            title(tla, sprintf('Cellpose sweep — %s', mTag), ...
                  'Color','w', 'FontSize',22, 'FontWeight','bold', 'Interpreter','none');
            ax1 = nexttile; heatmapPanel(ax1, cpVals, ftVals, nObj_mn, 'n objects',  'parula', '%d');
            ax2 = nexttile; heatmapPanel(ax2, cpVals, ftVals, cov_mn,  'coverage %', 'hot',    '%.1f');

            % Visual label grid
            figNum = figNum + 1;
            nCols = nFT + 1;
            tileW = min(280, floor(1760 / nCols));
            tileH = round(tileW * imgH / imgW);
            hfb = figure(figNum);
            set(hfb, 'Name', sprintf('Cellpose labels — %s', mTag), ...
                     'NumberTitle','off', 'Color','k', ...
                     'Position', [50 50 tileW*nCols+80 min(1050,tileH*nCP+80)]);
            tlb = tiledlayout(nCP, nCols, 'TileSpacing','none', 'Padding','tight');
            title(tlb, sprintf('Cellpose segmentation — %s', mTag), ...
                  'Color','w', 'FontSize',20, 'Interpreter','none');

            for ci = 1:nCP
                ax = nexttile; imshow(Ireal, []); hold on;
                text(4, 6, sprintf('cp=%+d', cpVals(ci)), ...
                     'Color','y', 'FontSize',14, 'FontWeight','bold', ...
                     'VerticalAlignment','top', 'Interpreter','none');
                if ci == 1, title('Raw','Color','w','FontSize',14,'FontWeight','normal'); end
                hold off;

                for fi = 1:nFT
                    ax = nexttile;
                    Lk = L_mn{ci,fi};
                    if max(Lk(:)) > 0
                        rgb = label2rgb(Lk, 'hsv', 'k', 'shuffle');
                    else
                        rgb = zeros(imgH, imgW, 3, 'uint8');
                    end
                    imshow(rgb); hold on;
                    n_k = nObj_mn(ci,fi); c_k = cov_mn(ci,fi);
                    if ~isnan(n_k)
                        text(4, 6, sprintf('n=%d  %.0f%%', n_k, c_k), ...
                             'Color','y', 'FontSize',12, ...
                             'VerticalAlignment','top', 'Interpreter','none');
                    end
                    if ci == 1
                        title(sprintf('ft=%.2f', ftVals(fi)), ...
                              'Color','w', 'FontSize',14, 'FontWeight','normal');
                    end
                    hold off;
                end
            end

            fprintf('  Figs %d-%d: heatmap + label grid for %s\n', figNum-1, figNum, mTag);

            if doExport && isfolder(outDir)
                safeMTag = strrep(strrep(strtrim(mTag), ' ', '_'), '=', '-');
                fn_heat = fullfile(outDir, sprintf('cellpose_sweep_heatmaps_%s.pdf', safeMTag));
                fn_labs = fullfile(outDir, sprintf('cellpose_sweep_labels_%s.pdf',  safeMTag));
                exportgraphics(hfa, fn_heat, 'ContentType', 'vector');
                exportgraphics(hfb, fn_labs, 'ContentType', 'vector');
                fn_tab = fullfile(outDir, sprintf('cellpose_sweep_table_%s.tex', safeMTag));
                writeCPtableLatex(fn_tab, cpModels{mi}, niStr, cpVals, ftVals, nObj_mn, cov_mn, cpMinSize);
                fprintf('  Exported: %s\n  Exported: %s\n  Exported: %s\n', fn_heat, fn_labs, fn_tab);
            end
        end
    end

    fprintf('  %d figures generated.\n', figNum);
end   % runCelloseSweep


% -------------------------------------------------------------------------
function runBlobFallback()
% Blob fallback: marker-controlled watershed sweep when Cellpose is absent.
% Sweeps threshold × hMinima; display mirrors the Cellpose figure style.

    nTH = numel(blobThreshVals);
    nHM = numel(blobHMinVals);

    % Build cache tag from image size + grid
    blobCacheTagNew = sprintf('blob|sz%dx%d|th%s|hm%s|ms%d', ...
        size(Ireal,1), size(Ireal,2), ...
        num2str(blobThreshVals,'%.2f '), ...
        num2str(blobHMinVals,'%.0f '), blobMinSize);

    doSweep = true;
    if evalin('base','exist(''blobCacheTag'',''var'')') && ...
       strcmp(evalin('base','blobCacheTag'), blobCacheTagNew) && ...
       evalin('base','exist(''blobLall'',''var'')')
        fprintf('=== Using cached blob results ===\n');
        fprintf('    Change a grid parameter or clear blobLall to force rerun.\n\n');
        doSweep    = false;
        blobLall    = evalin('base','blobLall');
        blobNObjAll = evalin('base','blobNObjAll');
        blobCovAll  = evalin('base','blobCovAll');
    end

    if doSweep
        blobNObjAll = nan(nTH, nHM);
        blobCovAll  = nan(nTH, nHM);
        blobLall    = cell(nTH, nHM);

        fprintf('=== Blob fallback sweep  (%d configurations) ===\n', nTH * nHM);

        for ti = 1:nTH
            for hi = 1:nHM
                fprintf('  thresh=%.2f  hMin=%d  ... ', blobThreshVals(ti), blobHMinVals(hi));
                tic;
                try
                    [~, Lk] = watershedSegment(Ireal, ...
                        'method',    'marker', ...
                        'threshold', blobThreshVals(ti), ...
                        'hMinima',   blobHMinVals(hi), ...
                        'minArea',   blobMinSize);
                    n_k = double(max(Lk(:)));
                    c_k = 100 * nnz(Lk > 0) / numel(Lk);
                    blobNObjAll(ti,hi) = n_k;
                    blobCovAll(ti,hi)  = c_k;
                    blobLall{ti,hi}    = Lk;
                    fprintf('%.1fs  n=%3d  cov=%5.1f%%\n', toc, n_k, c_k);
                catch ME
                    fprintf('FAILED (%s)\n', ME.message);
                    blobLall{ti,hi} = zeros(size(Ireal),'uint16');
                end
            end
        end

        assignin('base', 'blobLall',    blobLall);
        assignin('base', 'blobNObjAll', blobNObjAll);
        assignin('base', 'blobCovAll',  blobCovAll);
        assignin('base', 'blobCacheTag', blobCacheTagNew);
    end

    % ---- Console table --------------------------------------------------
    printBlobTable(blobThreshVals, blobHMinVals, blobNObjAll, blobCovAll);

    % ---- Heatmap figure -------------------------------------------------
    [imgH, imgW] = size(Ireal);
    hfa = figure(1);
    set(hfa, 'Name', 'Blob fallback heatmaps', ...
             'NumberTitle','off', 'Color',[0.13 0.13 0.13]);
    tla = tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','compact');
    title(tla, 'Blob fallback sweep (marker watershed)', ...
          'Color','w', 'FontSize',22, 'FontWeight','bold');
    ax1 = nexttile;
    heatmapPanel(ax1, blobThreshVals, blobHMinVals, blobNObjAll, ...
                 'n objects', 'parula', '%d');
    xlabel(ax1, 'hMinima',   'Color','w');
    ylabel(ax1, 'threshold', 'Color','w');
    ax2 = nexttile;
    heatmapPanel(ax2, blobThreshVals, blobHMinVals, blobCovAll, ...
                 'coverage %', 'hot', '%.1f');
    xlabel(ax2, 'hMinima',   'Color','w');
    ylabel(ax2, 'threshold', 'Color','w');

    % ---- Visual label grid ----------------------------------------------
    nCols = nHM + 1;
    tileW = min(280, floor(1760 / nCols));
    tileH = round(tileW * imgH / imgW);
    hfb = figure(2);
    set(hfb, 'Name', 'Blob fallback labels', ...
             'NumberTitle','off', 'Color','k', ...
             'Position', [50 50 tileW*nCols+80 min(1050,tileH*nTH+80)]);
    tlb = tiledlayout(nTH, nCols, 'TileSpacing','none', 'Padding','tight');
    title(tlb, 'Blob fallback — marker watershed', ...
          'Color','w', 'FontSize',20);

    for ti = 1:nTH
        ax = nexttile; imshow(Ireal, []); hold on;
        text(4, 6, sprintf('th=%.2f', blobThreshVals(ti)), ...
             'Color','y', 'FontSize',14, 'FontWeight','bold', ...
             'VerticalAlignment','top', 'Interpreter','none');
        if ti == 1, title('Raw','Color','w','FontSize',14,'FontWeight','normal'); end
        hold off;

        for hi = 1:nHM
            ax = nexttile;
            Lk = blobLall{ti,hi};
            if max(Lk(:)) > 0
                rgb = label2rgb(Lk, 'hsv', 'k', 'shuffle');
            else
                rgb = zeros(imgH, imgW, 3, 'uint8');
            end
            imshow(rgb); hold on;
            n_k = blobNObjAll(ti,hi); c_k = blobCovAll(ti,hi);
            if ~isnan(n_k)
                text(4, 6, sprintf('n=%d  %.0f%%', n_k, c_k), ...
                     'Color','y', 'FontSize',12, ...
                     'VerticalAlignment','top', 'Interpreter','none');
            end
            if ti == 1
                title(sprintf('hMin=%d', blobHMinVals(hi)), ...
                      'Color','w', 'FontSize',14, 'FontWeight','normal');
            end
            hold off;
        end
    end

    fprintf('  2 figures generated (heatmaps + label grid).\n');
end   % runBlobFallback


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function Iout = loadImage(wsVar, imgPath, imgFile)
% Load image: workspace variable takes priority over file.
    Iout = [];

    % 1. Try workspace variable
    if ~isempty(wsVar) && evalin('base', sprintf('exist(''%s'',''var'')', wsVar))
        raw = evalin('base', wsVar);
        if isnumeric(raw) && ~isempty(raw)
            if size(raw,3) > 1, raw = rgb2gray(raw); end
            Iout = im2single(squeeze(raw(:,:,1,1,1)));
            fprintf('=== Image from workspace variable ''%s''  (%d x %d px) ===\n\n', ...
                    wsVar, size(Iout,2), size(Iout,1));
            return;
        end
    end

    % 2. Fall back to file
    mitoMat = fullfile(imgPath, imgFile);
    if ~exist(mitoMat, 'file')
        if ~isempty(wsVar)
            error('sweepCellpose:noImage', ...
                  'Workspace variable ''%s'' not found and file not found:\n  %s', ...
                  wsVar, mitoMat);
        else
            error('sweepCellpose:noImage', 'Image file not found:\n  %s', mitoMat);
        end
    end
    tmp  = load(mitoMat);
    flds = fieldnames(tmp);
    raw  = tmp.(flds{1});
    if size(raw,3) > 1, raw = rgb2gray(raw); end
    Iout = im2single(raw);
    fprintf('=== Image loaded from file: %s  (%d x %d px,  var: %s) ===\n\n', ...
            imgFile, size(Iout,2), size(Iout,1), flds{1});
end


function printCPtable(model, niStr, cpVals, ftVals, nObjM, covM)
% Print 2-D table of n_objects / coverage% for one Cellpose model.
nCP  = numel(cpVals);
nFT  = numel(ftVals);
colW = 12;
sep  = repmat('=', 1, 10 + nFT * colW);
fprintf('%s\nCELLPOSE SWEEP — %s%s\n%s\n', sep, model, niStr, sep);
fprintf('  %-8s', 'cp \ ft');
for fi = 1:nFT, fprintf('  %*s', colW-2, sprintf('ft=%.2f', ftVals(fi))); end
fprintf('\n  %s\n', repmat('-', 1, 8 + nFT * colW));
for ci = 1:nCP
    fprintf('  cp=%+d  ', cpVals(ci));
    for fi = 1:nFT
        n_k = nObjM(ci,fi); c_k = covM(ci,fi);
        if isnan(n_k), fprintf('  %*s', colW-2, 'ERR');
        else,          fprintf('  %4d/%.1f%%', n_k, c_k); end
    end
    fprintf('\n');
end
fprintf('%s\n  (format: n_objects / coverage%%)\n\n', sep);
end


function printBlobTable(threshVals, hMinVals, nObjM, covM)
% Print 2-D table for the blob fallback sweep.
nTH  = numel(threshVals);
nHM  = numel(hMinVals);
colW = 12;
sep  = repmat('=', 1, 12 + nHM * colW);
fprintf('%s\nBLOB FALLBACK SWEEP (marker watershed)\n%s\n', sep, sep);
fprintf('  %-10s', 'th \ hMin');
for hi = 1:nHM, fprintf('  %*s', colW-2, sprintf('hMin=%d', hMinVals(hi))); end
fprintf('\n  %s\n', repmat('-', 1, 10 + nHM * colW));
for ti = 1:nTH
    fprintf('  th=%.2f  ', threshVals(ti));
    for hi = 1:nHM
        n_k = nObjM(ti,hi); c_k = covM(ti,hi);
        if isnan(n_k), fprintf('  %*s', colW-2, 'ERR');
        else,          fprintf('  %4d/%.1f%%', n_k, c_k); end
    end
    fprintf('\n');
end
fprintf('%s\n  (format: n_objects / coverage%%)\n\n', sep);
end


function heatmapPanel(ax, rowVals, colVals, data, titleStr, cmap, fmt)
% Draw an annotated heatmap: rows = rowVals, cols = colVals.
nR = numel(rowVals);
nC = numel(colVals);
imagesc(ax, 1:nC, 1:nR, data);
colormap(ax, cmap);
cb = colorbar(ax); cb.Color = 'w';
set(ax, ...
    'XTick', 1:nC, 'XTickLabel', arrayfun(@(v) sprintf('%.2g',v), colVals, 'UniformOutput',false), ...
    'YTick', 1:nR, 'YTickLabel', arrayfun(@(v) sprintf('%.2g',v), rowVals, 'UniformOutput',false), ...
    'Color', [0.1 0.1 0.1], 'XColor','w', 'YColor','w', ...
    'TickDir', 'out', 'FontSize',18);
title(ax, titleStr, 'Color','w', 'FontSize',20);
for ri = 1:nR
    for ci = 1:nC
        v = data(ri,ci);
        if ~isnan(v)
            text(ax, ci, ri, sprintf(fmt, v), ...
                 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                 'Color','w', 'FontSize',16, 'FontWeight','bold');
        end
    end
end
end


function writeCPtableLatex(fname, model, niStr, cpVals, ftVals, nObjM, covM, minSize)
% Write a LaTeX tabular of Cellpose sweep results to a .tex file.
nCP = numel(cpVals);
nFT = numel(ftVals);
fid = fopen(fname, 'w');
if fid < 0
    warning('sweepCellpose:tableWrite', 'Cannot write table: %s', fname);
    return;
end
fprintf(fid, '%% Auto-generated by sweepCellpose.m\n');
fprintf(fid, '%% Model: %s%s   minSize: %d px^2\n', model, strtrim(niStr), minSize);
fprintf(fid, '\\begin{tabular}{r|%s}\n\\hline\n', repmat('r', 1, nFT));
fprintf(fid, '\\texttt{cp} $\\backslash$ \\texttt{ft}');
for fi = 1:nFT, fprintf(fid, ' & \\texttt{%.2f}', ftVals(fi)); end
fprintf(fid, ' \\\\\n\\hline\n');
for ci = 1:nCP
    if cpVals(ci) > 0, fprintf(fid, '$+%d$', cpVals(ci));
    else,              fprintf(fid, '$%d$',  cpVals(ci)); end
    for fi = 1:nFT
        n_k = nObjM(ci,fi); c_k = covM(ci,fi);
        if isnan(n_k), fprintf(fid, ' & ---');
        else,          fprintf(fid, ' & $%d$ / $%.1f$\\%%', n_k, c_k); end
    end
    fprintf(fid, ' \\\\\n');
end
fprintf(fid, '\\hline\n\\end{tabular}\n');
fclose(fid);
end
