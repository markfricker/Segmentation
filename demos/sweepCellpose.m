%%
% sweepCellpose.m
%
% Sweep Cellpose hyperparameters to identify optimal settings for
% mitochondria segmentation ground truth.
%
% Unlike sweepSegmentation (which treats Cellpose as fixed GT and sweeps
% classical algorithms), this script is entirely about tuning Cellpose
% itself.  There is no external ground truth — the output metrics
% (object count, foreground coverage) and visual inspection are used
% to identify the 'Goldilocks' parameter set: not under-segmenting
% faint/small objects, but not over-merging or producing noise detections.
%
% STRATEGY
%   For each (model, nIter) combination:
%     Sweep CellThreshold (cellProb) × FlowErrorThreshold (flowThreshold).
%     Record n_objects and foreground coverage for each combo.
%   Output a 2-D console table per model, plus heatmap and visual-grid
%   figures for fast identification of the best parameter region.
%
% CACHING
%   Cellpose sweep results (Lall, nObjAll, covAll) are cached in the
%   workspace.  Re-running the script with the same sweep grid but
%   different display options skips the Cellpose calls entirely.
%   To force a full rerun, clear the workspace or change any grid
%   parameter (cpModels, cpVals, ftVals, niVals, cpDiam).
%
% METRICS
%   n objects   — number of distinct detected objects (= max label index).
%                 Too high → over-segmentation or noise.
%                 Too low  → missing objects or under-segmentation.
%   coverage %  — fraction of image pixels labelled as foreground (× 100).
%                 Should match the expected mitochondrial area fraction.
%
% FIGURES (per model × nIter combination)
%   Fig A — Combined heatmaps: n_objects (left) and coverage % (right).
%           Rows = cellProb values, Cols = flowThreshold values.
%   Fig B — Visual label grid (full image).
%           Row per cellProb, column per flowThreshold.
%           Annotated with n_objects and coverage%.
%
% CONSOLE
%   Per (model, nIter): 2-D grid table of n_objects / coverage%.
%
% CONFIGURATION
%   cpModels  — cell array of Cellpose model names to test.
%   cpDiam    — object diameter in pixels (fixed; tune separately).
%   cpVals    — CellThreshold (cellProb) sweep values.  Range: [-6, 6].
%               Lower = more sensitive; try [-4..1] for mitochondria.
%   ftVals    — FlowErrorThreshold sweep values.  Range: [0.1, 3].
%               Higher = more objects (rougher boundaries).
%   niVals    — nIter values.  0 = Cellpose default (~200).
%               2000 recommended for elongated structures (literature).
%
% REQUIREMENTS
%   cellposeEnhance.m on path (BlobFilters toolbox).
%   Medical Imaging Toolbox Interface for Cellpose Library Add-On.
%   Python environment with cellpose installed.

clc; close all;   % intentionally no 'clear' — workspace cache is preserved

% =========================================================================
% CONFIGURATION
% =========================================================================
blobFunctionPath = 'C:\Users\dops0035\Documents\Research\Matlab Projects\BlobFilters_sandbox\src';
blobImagePath    = 'C:\Users\dops0035\Documents\Research\Matlab Projects\BlobFilters_sandbox\demos';

% ---- Cellpose models and nIter to test ----------------------------------
% cyto3: general-purpose Cellpose model.
% bact_fluor_cp3: Cellpose 3 bacteria-fluorescence model; morphologically
%   similar to mitochondria so worth comparison.
% nIter=2000: literature recommendation for elongated structures (cyto3 only).
cpModels  = {'cyto3', 'bact_fluor_cp3'};
niVals    = [0, 2000];
cpDiam    = 10;           % expected object diameter in pixels (fixed)

% Combinations to exclude: each entry is {modelName, nIterValue}.
% bact_fluor_cp3+nIter=2000 is excluded — bact_fluor_cp3 already underperforms
% at default nIter so the high-nIter variant adds no useful information.
skipCombos = {{'bact_fluor_cp3', 2000}};

% ---- Sweep grid ---------------------------------------------------------
% Narrowed around sweet spot found on 2026-03-03 (cp=0, ft=0.8).
% For a broad initial survey restore: cpVals=[-4,-3,-2,-1,0,1,2], ftVals=[0.30,0.50,0.80,1.20]
cpVals    = [-1, 0, 1];          % CellThreshold values
ftVals    = [0.60, 0.80, 1.00]; % FlowErrorThreshold values
cpMinSize = 50;                  % px^2 — discard objects below this area;
                                 % equalises false-positive rate across models

% ---- Figure export -------------------------------------------------------
% Set doExport = true to save sweep figures as PDFs in the BlobFilters
% docs/figures directory for inclusion in BlobFilters_manual.tex.
% Leave false until outputs have been inspected and parameters are finalised.
doExport = false;
outDir   = fullfile(blobFunctionPath, '..', 'docs', 'figures');
% Note: docs/ lives at BlobFilters_sandbox/docs/, one level above src/.

% =========================================================================
% Path setup
% =========================================================================
demoDir = fileparts(mfilename('fullpath'));
addpath(fullfile(demoDir, '..', 'src'));
if ~isempty(blobFunctionPath) && exist(blobFunctionPath, 'dir')
    addpath(blobFunctionPath);
end
if ~exist('cellposeEnhance', 'file')
    error('sweepCellpose:noCellposeEnhance', ...
          'cellposeEnhance not found.  Check blobFunctionPath.');
end
if exist('cellpose', 'file') == 0
    error('sweepCellpose:noAddon', ...
          ['Cellpose add-on not installed.\n' ...
           'Install "Medical Imaging Toolbox Interface for Cellpose Library" ' ...
           'via Home > Add-Ons.']);
end

% =========================================================================
% Cache check — skip Cellpose calls if sweep results are already in workspace
% =========================================================================
nM  = numel(cpModels);
nNI = numel(niVals);
nCP = numel(cpVals);
nFT = numel(ftVals);

% Build a tag that uniquely identifies the current sweep grid
cpCacheTagNew = sprintf('%s|ni%s|d%.0f|cp%s|ft%s|ms%.0f', ...
    strjoin(cpModels, '_'), num2str(niVals,'%g '), cpDiam, ...
    num2str(cpVals,'%.2f '), num2str(ftVals,'%.2f '), cpMinSize);

doSweep = true;
if exist('cpCacheTag','var')  && strcmp(cpCacheTag, cpCacheTagNew) && ...
   exist('Lall','var')        && isequal(size(Lall), [nM nNI nCP nFT]) && ...
   exist('nObjAll','var')     && exist('Ireal','var')
    fprintf('=== Using cached sweep results ===\n');
    fprintf('    (grid unchanged — only figures will be regenerated)\n');
    fprintf('    To force a full rerun: clear Lall, or change a grid parameter.\n\n');
    doSweep = false;
end

% =========================================================================
% 1.  Load image  (skipped if cache valid)
% =========================================================================
if doSweep
    fprintf('=== Loading image ===\n');
    mitoMat = fullfile(blobImagePath, 'mitImagecrop.mat');

    if ~exist(mitoMat, 'file')
        error('sweepCellpose:noImage', ...
              'mitImagecrop.mat not found in:\n  %s\nCheck blobImagePath.', ...
              blobImagePath);
    end

    % Load defensively — accept any single 2-D numeric variable in the file
    tmp = load(mitoMat);
    flds = fieldnames(tmp);
    Iraw = tmp.(flds{1});
    if size(Iraw, 3) > 1, Iraw = rgb2gray(Iraw); end
    Ireal = im2single(Iraw);
    fprintf('  Loaded mitImagecrop.mat  (%d x %d px,  var: %s)\n', ...
            size(Ireal,2), size(Ireal,1), flds{1});

    % =====================================================================
    % 2.  Run sweep
    % =====================================================================
    nObjAll = nan(nM, nNI, nCP, nFT);
    covAll  = nan(nM, nNI, nCP, nFT);
    Lall    = cell(nM, nNI, nCP, nFT);

    nTotal = nM * nNI * nCP * nFT;
    fprintf('\n=== Cellpose parameter sweep  (%d configurations) ===\n', nTotal);

    for mi = 1:nM
        for ni = 1:nNI
            niStr = '';
            if niVals(ni) > 0, niStr = sprintf('  nIter=%d', niVals(ni)); end

            % Skip excluded (model, nIter) combinations
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
                        [~, Lk] = cellposeEnhance(Ireal, pCP);
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

    cpCacheTag = cpCacheTagNew;
end  % doSweep

% =========================================================================
% 3.  Console tables  (always printed)
% =========================================================================
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

% =========================================================================
% 4.  Figures: heatmaps + visual label grid  (always regenerated)
% =========================================================================
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

        nObj_mn = squeeze(nObjAll(mi,ni,:,:));   % [nCP × nFT]
        cov_mn  = squeeze(covAll(mi,ni,:,:));    % [nCP × nFT]
        L_mn    = squeeze(Lall(mi,ni,:,:));      % {nCP × nFT} cell

        % ---- Figure A: heatmaps (n_objects + coverage%) side-by-side ----
        figNum = figNum + 1;
        hfa = figure(figNum);
        set(hfa, 'Name', sprintf('Cellpose heatmaps — %s', mTag), ...
                 'NumberTitle','off', 'Color',[0.13 0.13 0.13]);

        tla = tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','compact');
        title(tla, sprintf('Cellpose sweep — %s', mTag), ...
              'Color','w', 'FontSize',11, 'FontWeight','bold');

        ax1 = nexttile;
        heatmapPanel(ax1, cpVals, ftVals, nObj_mn, 'n objects',  'parula', '%d');

        ax2 = nexttile;
        heatmapPanel(ax2, cpVals, ftVals, cov_mn,  'coverage %', 'hot',    '%.1f');

        % ---- Figure B: visual label grid (full image) -------------------
        figNum = figNum + 1;

        % Scale tile size to the image aspect ratio so panels aren't distorted
        nCols   = nFT + 1;          % raw column + one per flowThreshold
        tileW   = min(280, floor(1760 / nCols));
        tileH   = round(tileW * imgH / imgW);
        figW    = tileW * nCols + 80;
        figH    = min(1050, tileH * nCP  + 80);

        hfb = figure(figNum);
        set(hfb, 'Name', sprintf('Cellpose labels — %s', mTag), ...
                 'NumberTitle','off', 'Color','k', ...
                 'Position', [50 50 figW figH]);

        tlb = tiledlayout(nCP, nCols, 'TileSpacing','none', 'Padding','tight');
        title(tlb, sprintf('Labels (full image)  —  %s', mTag), ...
              'Color','w', 'FontSize',10);

        for ci = 1:nCP
            % Column 1: raw image annotated with cellProb value
            ax = nexttile; %#ok<NASGU>
            imshow(Ireal, []);
            hold on;
            text(4, 6, sprintf('cp=%+d', cpVals(ci)), ...
                 'Color','y', 'FontSize',7, 'FontWeight','bold', ...
                 'VerticalAlignment','top', 'Interpreter','none');
            if ci == 1
                title('Raw', 'Color','w', 'FontSize',7, 'FontWeight','normal');
            end
            hold off;

            % Columns 2…nFT+1: labelled results
            for fi = 1:nFT
                ax = nexttile; %#ok<NASGU>
                Lk = L_mn{ci,fi};
                if max(Lk(:)) > 0
                    rgb = label2rgb(Lk, 'hsv', 'k', 'shuffle');
                else
                    rgb = zeros(imgH, imgW, 3, 'uint8');
                end
                imshow(rgb);
                hold on;
                n_k = nObj_mn(ci,fi);
                c_k = cov_mn(ci,fi);
                if ~isnan(n_k)
                    text(4, 6, sprintf('n=%d  %.0f%%', n_k, c_k), ...
                         'Color','y', 'FontSize',6, ...
                         'VerticalAlignment','top', 'Interpreter','none');
                end
                if ci == 1
                    title(sprintf('ft=%.2f', ftVals(fi)), ...
                          'Color','w', 'FontSize',7, 'FontWeight','normal');
                end
                hold off;
            end
        end

        fprintf('  Figs %d-%d: heatmap + label grid for %s\n', ...
                figNum-1, figNum, mTag);

        % ---- Export figures for LaTeX manual ----------------------------
        if doExport && ~isempty(outDir) && isfolder(outDir)
            % Replace spaces and '=' so filenames are clean for LaTeX \includegraphics
            safeMTag = strrep(strrep(strtrim(mTag), ' ', '_'), '=', '-');
            fn_heat  = fullfile(outDir, ...
                sprintf('cellpose_sweep_heatmaps_%s.pdf', safeMTag));
            fn_labs  = fullfile(outDir, ...
                sprintf('cellpose_sweep_labels_%s.pdf', safeMTag));
            exportgraphics(hfa, fn_heat, 'ContentType', 'vector');
            exportgraphics(hfb, fn_labs, 'ContentType', 'vector');
            fprintf('  Exported: %s\n', fn_heat);
            fprintf('  Exported: %s\n', fn_labs);
        elseif doExport
            warning('sweepCellpose:noOutDir', ...
                    'doExport=true but outDir does not exist:\n  %s', outDir);
        end
    end
end

fprintf('\nDone.  %d figures generated.\n', figNum);


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function printCPtable(model, niStr, cpVals, ftVals, nObjM, covM)
% printCPtable  Print 2-D table of n_objects / coverage% for one model.
nCP = numel(cpVals);
nFT = numel(ftVals);
colW = 12;
sep  = repmat('=', 1, 10 + nFT * colW);

fprintf('%s\n', sep);
fprintf('CELLPOSE SWEEP — %s%s\n', model, niStr);
fprintf('%s\n', sep);
fprintf('  %-8s', 'cp \ ft');
for fi = 1:nFT
    fprintf('  %*s', colW-2, sprintf('ft=%.2f', ftVals(fi)));
end
fprintf('\n  %s\n', repmat('-', 1, 8 + nFT * colW));
for ci = 1:nCP
    fprintf('  cp=%+d  ', cpVals(ci));
    for fi = 1:nFT
        n_k = nObjM(ci,fi);
        c_k = covM(ci,fi);
        if isnan(n_k)
            fprintf('  %*s', colW-2, 'ERR');
        else
            fprintf('  %4d/%.1f%%', n_k, c_k);
        end
    end
    fprintf('\n');
end
fprintf('%s\n', sep);
fprintf('  (format: n_objects / coverage%%)\n\n');
end


function heatmapPanel(ax, cpVals, ftVals, data, titleStr, cmap, fmt)
% heatmapPanel  Draw an annotated heatmap: rows = cpVals, cols = ftVals.
nCP = numel(cpVals);
nFT = numel(ftVals);

imagesc(ax, 1:nFT, 1:nCP, data);
colormap(ax, cmap);
cb = colorbar(ax);
cb.Color = 'w';

set(ax, ...
    'XTick',      1:nFT, ...
    'XTickLabel', arrayfun(@(v) sprintf('%.2f',v), ftVals, 'UniformOutput',false), ...
    'YTick',      1:nCP, ...
    'YTickLabel', arrayfun(@(v) sprintf('%+d',v),  cpVals, 'UniformOutput',false), ...
    'Color',      [0.1 0.1 0.1], 'XColor','w', 'YColor','w', ...
    'TickDir',    'out', 'FontSize',9);
xlabel(ax, 'flowThreshold', 'Color','w');
ylabel(ax, 'cellProb',      'Color','w');
title(ax,  titleStr,        'Color','w', 'FontSize',10);

for ci = 1:nCP
    for fi = 1:nFT
        v = data(ci,fi);
        if ~isnan(v)
            text(ax, fi, ci, sprintf(fmt, v), ...
                 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                 'Color','w', 'FontSize',8, 'FontWeight','bold');
        end
    end
end
end
