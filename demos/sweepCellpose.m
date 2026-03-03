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
%   Fig B — Visual label grid (zoom region).
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
%   roi       — [x y w h] zoom crop for visual grid panels.
%
% REQUIREMENTS
%   cellposeEnhance.m on path (BlobFilters toolbox).
%   Medical Imaging Toolbox Interface for Cellpose Library Add-On.
%   Python environment with cellpose installed.

clear; clc; close all;

% =========================================================================
% CONFIGURATION
% =========================================================================
blobFunctionPath = 'C:\Users\dops0035\Documents\Research\Matlab Projects\BlobFilters_sandbox\src';
blobImagePath    = 'C:\Users\dops0035\Documents\Research\Matlab Projects\BlobFilters_sandbox\src\.claude\worktrees\thirsty-wescoff\demos';

% ---- Cellpose models and nIter to test ----------------------------------
% Add 'bact_fluor_cp3' to cpModels to compare models side-by-side.
% Add 2000 to niVals to compare default vs high-iteration mode.
cpModels  = {'cyto3'};    % e.g. {'cyto3', 'bact_fluor_cp3'}
niVals    = [0];          % 0 = Cellpose default; add 2000 to compare
cpDiam    = 10;           % expected object diameter in pixels (fixed)

% ---- Sweep grid ---------------------------------------------------------
% Broad initial sweep; narrow once the useful region is identified.
cpVals = [-4, -3, -2, -1, 0, 1, 2];   % CellThreshold values
ftVals = [0.30, 0.50, 0.80, 1.20];    % FlowErrorThreshold values

% ---- Display ------------------------------------------------------------
roi = [100  80  220 220];   % [x y w h] zoom region for visual grid panels

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
          'cellposeEnhance not found. Check blobFunctionPath.');
end
if exist('cellpose', 'file') == 0
    error('sweepCellpose:noAddon', ...
          ['Cellpose add-on not installed.\n' ...
           'Install "Medical Imaging Toolbox Interface for Cellpose Library" ' ...
           'via Home > Add-Ons.']);
end

% =========================================================================
% 1.  Load image
% =========================================================================
fprintf('=== Loading image ===\n');
mitoMat = fullfile(blobImagePath, 'mitImage.mat');
mitoImg = fullfile(blobImagePath, 'Ireal.png');

if exist(mitoMat, 'file')
    tmp   = load(mitoMat, 'I');
    Ireal = im2single(tmp.I);
    fprintf('  Loaded mitImage.mat  (%d x %d px)\n', size(Ireal,2), size(Ireal,1));
elseif exist(mitoImg, 'file')
    Ireal = im2single(imread(mitoImg));
    if size(Ireal,3) > 1, Ireal = rgb2gray(Ireal); end
    fprintf('  Loaded Ireal.png  (%d x %d px)\n', size(Ireal,2), size(Ireal,1));
else
    error('sweepCellpose:noImage', 'Image not found.  Check blobImagePath.');
end

% =========================================================================
% 2.  Run sweep
% =========================================================================
nM  = numel(cpModels);
nCP = numel(cpVals);
nFT = numel(ftVals);
nNI = numel(niVals);

nObjAll = nan(nM, nNI, nCP, nFT);   % object count per combo
covAll  = nan(nM, nNI, nCP, nFT);   % coverage % per combo
Lall    = cell(nM, nNI, nCP, nFT);  % label images (uint16)

nTotal = nM * nNI * nCP * nFT;
fprintf('\n=== Cellpose parameter sweep  (%d configurations) ===\n', nTotal);

for mi = 1:nM
    for ni = 1:nNI
        niStr = '';
        if niVals(ni) > 0, niStr = sprintf('  nIter=%d', niVals(ni)); end
        fprintf('\n--- %s%s ---\n', cpModels{mi}, niStr);
        t0 = tic;

        for ci = 1:nCP
            for fi = 1:nFT
                pCP.model         = cpModels{mi};
                pCP.diameter      = cpDiam;
                pCP.cellProb      = cpVals(ci);
                pCP.flowThreshold = ftVals(fi);
                pCP.nIter         = niVals(ni);

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

% =========================================================================
% 3.  Console tables
% =========================================================================
for mi = 1:nM
    for ni = 1:nNI
        niStr = '';
        if niVals(ni) > 0, niStr = sprintf('  nIter=%d', niVals(ni)); end

        fprintf('\n');
        printCPtable(cpModels{mi}, niStr, cpVals, ftVals, ...
                     squeeze(nObjAll(mi,ni,:,:)), ...
                     squeeze(covAll(mi,ni,:,:)));
    end
end

% =========================================================================
% 4.  Figures: heatmaps + visual label grids
% =========================================================================
figNum = 0;

for mi = 1:nM
    for ni = 1:nNI
        niStr = '';
        if niVals(ni) > 0, niStr = sprintf(' nIter=%d', niVals(ni)); end
        mTag = sprintf('%s%s', cpModels{mi}, niStr);

        nObj_mn = squeeze(nObjAll(mi,ni,:,:));  % [nCP × nFT]
        cov_mn  = squeeze(covAll(mi,ni,:,:));   % [nCP × nFT]
        L_mn    = squeeze(Lall(mi,ni,:,:));     % {nCP × nFT} cell

        % ---- Figure A: heatmaps side-by-side ----------------------------
        figNum = figNum + 1;
        hfa = figure(figNum);
        set(hfa, 'Name', sprintf('Cellpose heatmaps — %s', mTag), ...
                 'NumberTitle','off', 'Color', [0.13 0.13 0.13]);

        tla = tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
        title(tla, sprintf('Cellpose sweep — %s', mTag), ...
              'Color','w', 'FontSize',11, 'FontWeight','bold');

        % Left: n_objects
        ax1 = nexttile;
        heatmapPanel(ax1, cpVals, ftVals, nObj_mn, 'n objects', 'parula', '%d');

        % Right: coverage %
        ax2 = nexttile;
        heatmapPanel(ax2, cpVals, ftVals, cov_mn, 'coverage %', 'hot', '%.1f');

        % ---- Figure B: visual label grid ---------------------------------
        figNum = figNum + 1;

        % Dynamic figure size: ~170 px per tile
        tilePx  = 170;
        nCols   = nFT + 1;   % 1 raw column + nFT label columns
        figW    = min(1800, tilePx * nCols + 80);
        figH    = min(1100, tilePx * nCP  + 80);

        hfb = figure(figNum);
        set(hfb, 'Name', sprintf('Cellpose labels — %s', mTag), ...
                 'NumberTitle','off', 'Color','k', ...
                 'Position', [70 50 figW figH]);

        tlb = tiledlayout(nCP, nCols, 'TileSpacing','none','Padding','tight');
        title(tlb, sprintf('Labels (zoom)  —  %s', mTag), ...
              'Color','w', 'FontSize',10);

        Izoom = imcrop(Ireal, roi);

        for ci = 1:nCP
            % Column 1: raw image with cp label
            ax = nexttile;
            imshow(Izoom, []);
            hold on;
            text(ax, 4, 6, sprintf('cp=%+d', cpVals(ci)), ...
                 'Color','y', 'FontSize',7, 'FontWeight','bold', ...
                 'VerticalAlignment','top', 'Interpreter','none');
            if ci == 1
                title('Raw','Color','w','FontSize',7,'FontWeight','normal');
            end
            hold off;

            % Columns 2..nFT+1: labeled results
            for fi = 1:nFT
                ax = nexttile;
                Lk = L_mn{ci,fi};
                Lk_z = imcrop(Lk, roi);
                if max(Lk_z(:)) > 0
                    rgb = label2rgb(Lk_z, 'hsv', 'k', 'shuffle');
                else
                    rgb = zeros([size(Lk_z,1), size(Lk_z,2), 3], 'uint8');
                end
                imshow(rgb);
                hold on;
                n_k = nObj_mn(ci,fi);
                c_k = cov_mn(ci,fi);
                if ~isnan(n_k)
                    text(ax, 4, 6, sprintf('n=%d  %.0f%%', n_k, c_k), ...
                         'Color','y', 'FontSize',6, ...
                         'VerticalAlignment','top', 'Interpreter','none');
                end
                if ci == 1
                    title(sprintf('ft=%.2f', ftVals(fi)), ...
                          'Color','w','FontSize',7,'FontWeight','normal');
                end
                hold off;
            end
        end

        fprintf('  Figs %d-%d: heatmap + label grid for %s\n', ...
                figNum-1, figNum, mTag);
    end
end

fprintf('\nDone.  %d figures generated.\n', figNum);


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function printCPtable(model, niStr, cpVals, ftVals, nObjM, covM)
% printCPtable  Print 2-D table of n_objects / coverage% for one model.
%   Rows = cpVals, Cols = ftVals.  Format: nObj/cov%
nCP = numel(cpVals);
nFT = numel(ftVals);

sep = repmat('=', 1, 12 + nFT*12);
fprintf('%s\n', sep);
fprintf('CELLPOSE SWEEP — %s%s\n', model, niStr);
fprintf('%s\n', sep);
fprintf('  %-8s', 'cp \ ft');
for fi = 1:nFT
    fprintf('  %10s', sprintf('ft=%.2f', ftVals(fi)));
end
fprintf('\n  %s\n', repmat('-', 1, 10 + nFT*12));
for ci = 1:nCP
    fprintf('  cp=%+d  ', cpVals(ci));
    for fi = 1:nFT
        n_k = nObjM(ci,fi);
        c_k = covM(ci,fi);
        if isnan(n_k)
            fprintf('  %10s', 'ERR');
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
% heatmapPanel  Draw an annotated heatmap on axes ax.
%   Rows = cpVals, Cols = ftVals.
nCP = numel(cpVals);
nFT = numel(ftVals);

imagesc(ax, 1:nFT, 1:nCP, data);
colormap(ax, cmap);
cb = colorbar(ax);
cb.Color = 'w';

set(ax, 'XTick', 1:nFT, ...
        'XTickLabel', arrayfun(@(v)sprintf('%.2f',v), ftVals, 'UniformOutput',false), ...
        'YTick', 1:nCP, ...
        'YTickLabel', arrayfun(@(v)sprintf('%+d',v), cpVals, 'UniformOutput',false), ...
        'Color', [0.1 0.1 0.1], 'XColor','w', 'YColor','w', ...
        'GridColor','w', 'TickDir','out', 'FontSize',9);
xlabel(ax, 'flowThreshold', 'Color','w');
ylabel(ax, 'cellProb',      'Color','w');
title(ax, titleStr, 'Color','w', 'FontSize',10);

% Annotate each cell with its value
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
