%%
% demoSegmentation.m
%
% Demo script: enhancer parameter sweep on a real mitochondrial image.
%
% PIPELINE
%   Load mitImagecrop.mat from BlobFilters demos/.
%   Preprocessing : orientedGaussSmooth  → Ism
%   Cellpose GT   : cyto3, cellProb=0, flowThreshold=0.8, minSize=64
%   For each enhancer (log, fiber, capS, capD, rod):
%     Sweep enhancer parameter configs x watershed/hysteresis combos.
%     One figure per enhancer: Raw | best-WS per cfg | Cellpose GT.
%
% FIGURES PRODUCED
%   Fig 1 — logEnhance sweep        (best WS per config vs GT)
%   Fig 2 — fiberEnhance sweep
%   Fig 3 — capsule(single) sweep
%   Fig 4 — capsule(DoC) sweep      (3-row grid: wideWidth x alpha)
%   Fig 5 — rodGranulometry sweep
%   Fig 6 — HSV overlay summary     (hue=label, value=image intensity)
%
% REQUIREMENTS
%   localThresholdFast.m and watershedSegment.m (Segmentation_sandbox/src)
%   BlobFilters: logEnhance, fiberEnhance, capsuleEnhance,
%                rodGranulometryEnhance, cellposeEnhance, orientedGaussSmooth
%   Image Processing Toolbox (bwdist, watershed, imextendedmax, etc.)
%   MATLAB Medical Imaging Toolbox Interface for Cellpose Library

clear; clc; close all;

% =========================================================================
% CONFIGURATION  ← adjust these paths/flags for your machine
% =========================================================================
% Directory containing BlobFilters enhancer .m files.
blobFunctionPath = 'C:\Users\dops0035\Documents\Research\Matlab Projects\BlobFilters_sandbox\src';

% Set true to export all figures as PDF to outDir (for LaTeX manual).
doExport = false;
outDir   = fullfile(blobFunctionPath, '..', 'docs', 'figures');

% =========================================================================
% Path setup
% =========================================================================
demoDir = fileparts(mfilename('fullpath'));
addpath(fullfile(demoDir, '..', 'src'));

if ~isempty(blobFunctionPath) && exist(blobFunctionPath, 'dir')
    addpath(blobFunctionPath);
    hasBlobFilters = exist('rodGranulometryEnhance', 'file') ~= 0;
else
    hasBlobFilters = false;
end

if ~hasBlobFilters
    fprintf('BlobFilters not found at:\n  %s\nCheck blobFunctionPath.\n', blobFunctionPath);
    return
end

% =========================================================================
% 1.  Load image
% =========================================================================
fprintf('=== Loading image ===\n');
mitoMat = fullfile(blobFunctionPath, '..', 'demos', 'mitImagecrop.mat');
if ~exist(mitoMat, 'file')
    fprintf('mitImagecrop.mat not found at:\n  %s\n', mitoMat);
    return
end
tmp   = load(mitoMat);
flds  = fieldnames(tmp);
Ireal = im2single(tmp.(flds{1}));
fprintf('  Size: %d x %d px\n', size(Ireal,2), size(Ireal,1));

% =========================================================================
% 2.  Cellpose ground truth
% =========================================================================
fprintf('\n=== Cellpose ground truth ===\n');
cpGTmodel         = 'cyto3';
cpGTdiameter      = 10;
cpGTcellProb      = 0;
cpGTflowThreshold = 0.8;
cpGTnIter         = 0;
cpMinSize         = 64;

pCP.model         = cpGTmodel;
pCP.diameter      = cpGTdiameter;
pCP.cellProb      = cpGTcellProb;
pCP.flowThreshold = cpGTflowThreshold;
pCP.nIter         = cpGTnIter;
pCP.minSize       = cpMinSize;

if exist('cellpose','file') == 0
    error('demoSegmentation:noCellpose', ...
          'Cellpose add-on not installed. Cannot generate ground truth.');
end

fprintf('  cellposeEnhance (%s, cp=%d, ft=%.1f, minSize=%d)... ', ...
        cpGTmodel, cpGTcellProb, cpGTflowThreshold, cpMinSize);
tic;
[~, cpL] = cellposeEnhance(Ireal, pCP);
fprintf('%.1fs  (%d objects)\n', toc, max(cpL(:)));

cpBW  = logical(cpL > 0);
cpRGB = labelRGB(cpL);
nCP   = max(cpL(:));
gtLabel = sprintf('%s cp=%d ft=%.1f minSize=%d', ...
                  cpGTmodel, cpGTcellProb, cpGTflowThreshold, cpMinSize);

% =========================================================================
% 3.  OGS preprocessing
% =========================================================================
fprintf('\n=== OGS preprocessing ===\n');
pOGS.sigmaAlong = 4;  pOGS.sigmaAcross = 1.5;
pOGS.orientations = 8;  pOGS.sigmaGrad = 1.5;  pOGS.sigmaInt = 5;
fprintf('  orientedGaussSmooth... '); tic;
Ism = orientedGaussSmooth(Ireal, pOGS);
fprintf('%.1fs\n', toc);

% =========================================================================
% 4.  Watershed sweep grids (per-task threshold ranges)
% =========================================================================
hMinimaVals = [0.02 0.05 0.10];
minAreaVals = [20 50];

% Default: log, fiber, rod
wsThreshDef     = [0.20 0.25 0.30 0.35 0.40];
wsThreshLowDef  = {[], 0.10, 0.15};

% Capsule: higher thresholds, raised hysteresis lower bounds
wsThreshHigh    = [0.30 0.35 0.40 0.45 0.50];
wsThreshLowHigh = {[], 0.20, 0.25};

% =========================================================================
% 5.  Enhancer task definitions + per-task WS config
% =========================================================================
tasks  = buildTasks();
nTasks = numel(tasks);

for ti = 1:nTasks
    switch tasks{ti}.tag
        case {'capS-WS', 'capD-WS'}
            tasks{ti}.wsThresholds = wsThreshHigh;
            tasks{ti}.wsThreshLow  = wsThreshLowHigh;
        otherwise
            tasks{ti}.wsThresholds = wsThreshDef;
            tasks{ti}.wsThreshLow  = wsThreshLowDef;
    end
end

% =========================================================================
% 6.  Per-enhancer sweep — one figure per enhancer (Figs 1..nTasks)
% =========================================================================
fprintf('\n=== Enhancer parameter sweep ===\n');

overallBestLabels = cell(1, nTasks);   % best label image per enhancer (for HSV fig)
overallBestF1     = zeros(1, nTasks);  % corresponding F1

for ti = 1:nTasks
    tk = tasks{ti};

    % Build task-specific WS grid
    wsGridTask = {};
    for thr = tk.wsThresholds
        for hm = hMinimaVals
            for ma = minAreaVals
                for tli = 1:numel(tk.wsThreshLow)
                    tl = tk.wsThreshLow{tli};
                    if ~isempty(tl) && tl >= thr, continue; end
                    wsGridTask{end+1} = struct('threshold',thr,'hMinima',hm, ...
                                               'minArea',ma,'threshLow',tl); %#ok<SAGROW>
                end
            end
        end
    end
    nWSTask = numel(wsGridTask);
    nCfgs   = numel(tk.cfgs);

    fprintf('\n--- %s (%d configs x %d WS = %d evals) ---\n', ...
            tk.name, nCfgs, nWSTask, nCfgs * nWSTask);

    bestLabelsPerCfg = cell(1, nCfgs);
    bestF1PerCfg     = zeros(1, nCfgs);
    bestDescPerCfg   = cell(1, nCfgs);

    for ci = 1:nCfgs
        fprintf('  Config %d/%d: %s\n    Enhancement... ', ci, nCfgs, tk.cfgDescs{ci});
        tic;
        try
            R = tk.fn(Ism, tk.cfgs{ci});
            fprintf('%.1fs\n', toc);
        catch ME
            fprintf('SKIPPED (%s)\n', ME.message);
            bestLabelsPerCfg{ci} = zeros(size(Ireal), 'uint16');
            bestDescPerCfg{ci}   = sprintf('%s\n(SKIPPED)', tk.cfgDescs{ci});
            continue;
        end

        bestF1 = 0;
        bestL  = zeros(size(Ireal), 'uint16');
        bestWD = '';

        for wi = 1:nWSTask
            ws = wsGridTask{wi};
            try
                [BW, L] = watershedSegment(R, ...
                    'method',      'marker', ...
                    'threshold',   ws.threshold, ...
                    'threshLow',   ws.threshLow, ...
                    'smoothSigma', 1.5, ...
                    'hMinima',     ws.hMinima, ...
                    'minArea',     ws.minArea);
            catch
                continue;
            end
            [f1, ~, ~, ~] = evalBinary(BW, cpBW);
            if f1 > bestF1
                bestF1 = f1;
                bestL  = L;
                bestWD = wsDesc(ws);
            end
        end

        bestLabelsPerCfg{ci} = bestL;
        bestF1PerCfg(ci)     = bestF1;
        bestDescPerCfg{ci}   = sprintf('%s\nF1=%.3f  n=%d', ...
                                       tk.cfgDescs{ci}, bestF1, max(bestL(:)));
        fprintf('    Best F1=%.4f  WS: %s\n', bestF1, bestWD);
    end

    [bestF1overall, iBest] = max(bestF1PerCfg);

    overallBestLabels{ti} = bestLabelsPerCfg{iBest};
    overallBestF1(ti)     = bestF1overall;

    figNum   = ti;
    figTitle = sprintf('%s — best WS per config  |  best: %s, F1=%.3f  |  GT: %s', ...
                       tk.name, tk.cfgDescs{iBest}, bestF1overall, gtLabel);

    if strcmp(tk.tag, 'capD-WS')
        % DoC sweep: block layout (rows=wideWidth, cols=Raw+alphas+GT)
        docGridFigure(figNum, Ireal, cpRGB, bestLabelsPerCfg, bestF1PerCfg, ...
                      bestDescPerCfg, tk, gtLabel, nCP);
    else
        % Standard layout: Raw | cfg1 | cfg2 | ... | GT  (single row)
        panels      = [{Ireal}, cellfun(@labelRGB, bestLabelsPerCfg, 'UniformOutput',false), {cpRGB}];
        panelTitles = [{'Raw'}, bestDescPerCfg, {sprintf('Cellpose GT\n(n=%d)', nCP)}];
        cmaps       = [{'gray'}, repmat({[]}, 1, nCfgs), {[]}];
        figW        = max(900, 250 * (nCfgs + 2));
        figure(figNum);
        set(gcf, 'Name', sprintf('%s sweep', tk.name), 'NumberTitle','off', ...
                 'Color','k', 'Position',[30, 50, figW, 360]);
        segFillFigure(figNum, panels, figTitle, panelTitles, cmaps);
    end

    if doExport
        fname = sprintf('demoSeg_%s.pdf', strrep(tk.tag, '-', '_'));
        exportgraphics(figure(figNum), fullfile(outDir, fname), 'ContentType','vector');
        fprintf('  Exported: %s\n', fname);
    end

    fprintf('  >> Best overall: config %d (%s)  F1=%.4f\n', ...
            iBest, tk.cfgDescs{iBest}, bestF1overall);
    fprintf('  Fig %d generated.\n', figNum);
end

% =========================================================================
% HSV overlay summary figure  (Fig nTasks+1)
%   Hue = per-label colour, Saturation = 1, Value = original image intensity.
%   Background (unlabelled) shown as grayscale (S=0).
% =========================================================================
hsvPanels = [{labelHSVoverlay(zeros(size(Ireal),'uint16'), Ireal)}, ...   % raw (all grey)
             cellfun(@(L) labelHSVoverlay(L, Ireal), overallBestLabels, ...
                     'UniformOutput', false), ...
             {labelHSVoverlay(cpL, Ireal)}];
hsvTitles = [{'Raw'}];
for ti = 1:nTasks
    hsvTitles{end+1} = sprintf('%s\nF1=%.3f', tasks{ti}.tag, overallBestF1(ti)); %#ok<SAGROW>
end
hsvTitles{end+1} = sprintf('Cellpose GT\n(n=%d)', nCP);

figHSV = nTasks + 1;
figure(figHSV);
set(gcf, 'Name','HSV overlay — hue=label, value=image', 'NumberTitle','off', ...
         'Color','k', 'Position',[30, 50, 300*(nTasks+2), 360]);
segFillFigure(figHSV, hsvPanels, ...
    sprintf('HSV overlay: hue=label  value=image intensity  |  GT: %s', gtLabel), ...
    hsvTitles, repmat({[]}, 1, nTasks+2));

if doExport
    exportgraphics(figure(figHSV), fullfile(outDir, 'demoSeg_HSV_overlay.pdf'), ...
                   'ContentType','vector');
    fprintf('  Exported: demoSeg_HSV_overlay.pdf\n');
end

fprintf('Fig %d (HSV overlay) generated.\n', figHSV);
fprintf('\nDone. %d figures generated.\n', nTasks + 1);


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function rgb = labelRGB(L)
if max(L(:)) == 0
    rgb = zeros(size(L,1), size(L,2), 3, 'uint8');
else
    rgb = label2rgb(L, 'jet', 'k', 'shuffle');
end
end


function rgb = labelHSVoverlay(L, I)
% labelHSVoverlay  Colour-coded segmentation with original image intensity.
%   H = per-label hue (golden-ratio sequence, fixed seed for reproducibility)
%   S = 1 for labelled pixels, 0 for background (shown as grayscale)
%   V = normalised image intensity

% Normalise image to [0,1] for the Value channel
Iv = double(I);
Iv = (Iv - min(Iv(:))) / max(eps, max(Iv(:)) - min(Iv(:)));

nLabels = double(max(L(:)));
H = zeros(size(L), 'double');
S = zeros(size(L), 'double');   % background stays S=0 (grey)

if nLabels > 0
    % Golden-ratio hue sequence gives good perceptual separation
    rng(42);
    hues = mod((1:nLabels) * 0.618033988749895, 1);
    hues = hues(randperm(nLabels));

    % Vectorised assignment via lookup table
    hueMap = [0;     hues(:)];          % index 1 = label 0 (background)
    satMap = [0; ones(nLabels, 1)];
    idx    = double(L(:)) + 1;          % label 0 -> index 1
    H      = reshape(hueMap(idx), size(L));
    S      = reshape(satMap(idx), size(L));
end

rgb = im2uint8(hsv2rgb(cat(3, H, S, Iv)));
end


function segFillFigure(figNum, panels, figTitle, panelTitles, cmaps)
figure(figNum);
N = numel(panels);
for k = 1:N
    ax = subplot(1, N, k);
    if isempty(cmaps{k})
        imshow(panels{k});
    else
        imshow(panels{k}, []);
        colormap(ax, cmaps{k});
    end
    title(panelTitles{k}, 'Color','w', 'FontSize',9, 'Interpreter','none');
    set(ax, 'XColor','none', 'YColor','none');
end
sgtitle(figTitle, 'Color','w', 'FontSize',11, 'FontWeight','bold', 'Interpreter','none');
end


function docGridFigure(figNum, Ireal, cpRGB, bestLabelsPerCfg, bestF1PerCfg, ...
                       bestDescPerCfg, tk, gtLabel, nCP)
% docGridFigure  Block figure for 2-D parameter sweeps (e.g. DoC wideWidth x alpha).
%   Rows = tk.gridRows values (e.g. wideWidth), columns = Raw + tk.gridCols + GT.
nWW   = numel(tk.gridRows);
nAlph = numel(tk.gridCols);
nCols = nAlph + 2;             % Raw | alpha_1 .. alpha_n | GT

[~, iBest] = max(bestF1PerCfg);

figure(figNum);
set(gcf, 'Name', sprintf('%s sweep', tk.name), 'NumberTitle','off', ...
         'Color','k', 'Position',[30, 50, 260*nCols, 260*nWW + 60]);

tl = tiledlayout(nWW, nCols, 'TileSpacing','tight', 'Padding','compact');
title(tl, sprintf('%s — wideWidth (rows) x alpha (cols)  |  best: %s, F1=%.3f  |  GT: %s', ...
                  tk.name, tk.cfgDescs{iBest}, bestF1PerCfg(iBest), gtLabel), ...
     'Color','w', 'FontSize',11, 'FontWeight','bold', 'Interpreter','none');

for wi = 1:nWW
    % Column 1: Raw image (labelled with wideWidth for this row)
    ax = nexttile;
    imshow(Ireal, []);  colormap(ax, 'gray');
    title(ax, sprintf('Raw\nwW=%d', tk.gridRows(wi)), ...
          'Color','w', 'FontSize',9, 'Interpreter','none');
    set(ax, 'XColor','none', 'YColor','none');

    % Columns 2..nAlph+1: one per alpha value
    for ai = 1:nAlph
        ci  = (wi-1)*nAlph + ai;   % config index (wW outer, alpha inner)
        ax  = nexttile;
        imshow(labelRGB(bestLabelsPerCfg{ci}));
        if wi == 1
            ttl = sprintf('a=%.1f\nF1=%.3f  n=%d', ...
                          tk.gridCols(ai), bestF1PerCfg(ci), max(bestLabelsPerCfg{ci}(:)));
        else
            ttl = sprintf('F1=%.3f  n=%d', ...
                          bestF1PerCfg(ci), max(bestLabelsPerCfg{ci}(:)));
        end
        title(ax, ttl, 'Color','w', 'FontSize',9, 'Interpreter','none');
        set(ax, 'XColor','none', 'YColor','none');
    end

    % Last column: Cellpose GT
    ax = nexttile;
    imshow(cpRGB);
    if wi == 1
        title(ax, sprintf('Cellpose GT\n(n=%d)', nCP), ...
              'Color','w', 'FontSize',9, 'Interpreter','none');
    else
        title(ax, sprintf('(n=%d)', nCP), 'Color','w', 'FontSize',9, 'Interpreter','none');
    end
    set(ax, 'XColor','none', 'YColor','none');
end
end


function tasks = buildTasks()
% Configs tuned for mitochondria: width ~10 px, length ~15-50 px.
tasks = {};

% --- logEnhance -----------------------------------------------------------
% LoG optimal sigma = radius/sqrt(2); cross-section radius ~5 px -> sigma ~3-4.
t.name = 'logEnhance';  t.tag = 'log-WS';
t.fn   = @(I, p) logEnhance(I, p);
t.cfgs = {
    struct('sigmas',[2 3 4],   'normalize',true),
    struct('sigmas',[3 4 5],   'normalize',true),
    struct('sigmas',[4 5 6],   'normalize',true),
    struct('sigmas',[3 4 5 6], 'normalize',true)
};
t.cfgDescs = {
    'sigmas=[2 3 4]',
    'sigmas=[3 4 5]',
    'sigmas=[4 5 6]',
    'sigmas=[3 4 5 6]'
};
t.maxExpand = 4;
tasks{end+1} = t;

% --- fiberEnhance ---------------------------------------------------------
% widths = short-axis width; target ~10 px.
t.name = 'fiberEnhance';  t.tag = 'fib-WS';
t.fn   = @(I, p) fiberEnhance(I, p);
t.cfgs = {
    struct('widths',[8 9 10],       'multimode','stack','normalize',true),
    struct('widths',[9 10 11 12],   'multimode','stack','normalize',true),
    struct('widths',[10 12 14],     'multimode','stack','normalize',true)
};
t.cfgDescs = {
    'widths=[8 9 10]',
    'widths=[9 10 11 12]',
    'widths=[10 12 14]'
};
t.maxExpand = 7;
tasks{end+1} = t;

% --- capsuleEnhance single ------------------------------------------------
% lengths = half-lengths; full length 15-50 px -> half-lengths 8-25 px.
% orientations=12 (15 deg step) reduces angular discretisation artefacts.
t.name = 'capsule(single)';  t.tag = 'capS-WS';
t.fn   = @(I, p) capsuleEnhance(I, p);
t.cfgs = {
    struct('lengths',[8 12 16],    'width', 8,'orientations',12,'mode','single','normalize',true),
    struct('lengths',[12 16 20],   'width',10,'orientations',12,'mode','single','normalize',true),
    struct('lengths',[16 20 24],   'width',10,'orientations',12,'mode','single','normalize',true)
};
t.cfgDescs = {
    'lens=[8 12 16] w=8 or=12',
    'lens=[12 16 20] w=10 or=12',
    'lens=[16 20 24] w=10 or=12'
};
t.maxExpand = 5;
tasks{end+1} = t;

% --- capsuleEnhance DoC — sweep wideWidth x alpha -------------------------
% Fixed: width=10, lengths=[12 16 20 24], orientations=12.
% Sweep: wideWidth in [16 20 24] (1.6x, 2.0x, 2.4x inner width)
%        alpha     in [0.4 0.5 0.6 0.7]  (surround suppression weight)
% -> 12 combos total.
t.name = 'capsule(DoC)';  t.tag = 'capD-WS';
t.fn   = @(I, p) capsuleEnhance(I, p);
wideWidths = [16 20 24];
alphas     = [0.4 0.5 0.6 0.7];
t.cfgs     = {};
t.cfgDescs = {};
for ww = wideWidths
    for al = alphas
        t.cfgs{end+1} = struct('lengths',[12 16 20 24],'width',10, ...
            'wideWidth',ww,'alpha',al,'orientations',12,'mode','doc','normalize',true);
        t.cfgDescs{end+1} = sprintf('wW=%d a=%.1f', ww, al);
    end
end
t.gridRows  = wideWidths;   % rows of block figure (one per wideWidth value)
t.gridCols  = alphas;       % cols of block figure (one per alpha value)
t.maxExpand = 5;
tasks{end+1} = t;

% --- rodGranulometryEnhance -----------------------------------------------
% lengths = full line lengths; span measured range 15-50 px.
% orientations=12 (15 deg step).
t = struct();   % reset — prevents DoC-specific fields (gridRows/gridCols) carrying over
t.name = 'rodGranulometry';  t.tag = 'rod-WS';
t.fn   = @(I, p) rodGranulometryEnhance(I, p);
t.cfgs = {
    struct('lengths',[16 20 28],        'orientations',12,'normalize',true),
    struct('lengths',[20 28 36 48],     'orientations',12,'normalize',true),
    struct('lengths',[16 20 28 36 48],  'orientations',12,'normalize',true)
};
t.cfgDescs = {
    'lens=[16 20 28] or=12',
    'lens=[20 28 36 48] or=12',
    'lens=[16 20 28 36 48] or=12'
};
t.maxExpand = 8;
tasks{end+1} = t;
end


function [f1, prec, rec, iou] = evalBinary(BW, GT)
BW = logical(BW);  GT = logical(GT);
TP = nnz( BW &  GT);
FP = nnz( BW & ~GT);
FN = nnz(~BW &  GT);
prec = TP / max(1, TP + FP);
rec  = TP / max(1, TP + FN);
f1   = 2 * prec * rec / max(eps, prec + rec);
iou  = TP / max(1, TP + FP + FN);
end


function s = wsDesc(ws)
if isempty(ws.threshLow)
    s = sprintf('thr=%.2f h=%.2f minA=%d', ws.threshold, ws.hMinima, ws.minArea);
else
    s = sprintf('thr=%.2f lo=%.2f h=%.2f minA=%d', ...
                ws.threshold, ws.threshLow, ws.hMinima, ws.minArea);
end
end
