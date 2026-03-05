%%
% demoSegmentation.m
%
% Demo script: comparison of all segmentation methods on a synthetic test
% image and on a real mitochondrial fluorescence image.
%
% PIPELINE — Synthetic (Figs 1-3)
%   makeSyntheticBlobs()  — two touching blobs, elongated rod, small dot.
%   localThresholdFast    — sauvola, bernsen, bradley, niblack
%   watershedSegment      — distance, gradient, marker
%
% PIPELINE — Real mitochondrial image (Figs 4-8, one per enhancer)
%   Load mitImagecrop.mat from BlobFilters demos/.
%   Preprocessing : orientedGaussSmooth  → Ism
%   Cellpose GT   : cyto3, cellProb=0, flowThreshold=0.8, minSize=64
%   For each enhancer (log, fiber, capS, capD, rod):
%     Sweep enhancer parameter configs x watershed/hysteresis combos.
%     One figure per enhancer: Raw | best-WS per cfg | Cellpose GT.
%
% FIGURES PRODUCED
%   Fig 1 — Synthetic: thresholding  (Raw | Sauvola | Bernsen | Bradley | Niblack)
%   Fig 2 — Synthetic: watershed     (Raw | Sauvola mask | Dist | Grad | Marker)
%   Fig 3 — Synthetic: zoom on marker-WS
%   Fig 4 — logEnhance sweep        (best WS per config vs GT)
%   Fig 5 — fiberEnhance sweep
%   Fig 6 — capsule(single) sweep
%   Fig 7 — capsule(DoC) sweep
%   Fig 8 — rodGranulometry sweep
%
% REQUIREMENTS
%   localThresholdFast.m and watershedSegment.m (Segmentation_sandbox/src)
%   BlobFilters: logEnhance, fiberEnhance, capsuleEnhance,
%                rodGranulometryEnhance, cellposeEnhance, orientedGaussSmooth
%   Image Processing Toolbox (bwdist, watershed, imextendedmax, etc.)
%   MATLAB Medical Imaging Toolbox Interface for Cellpose Library

clear; clc; close all;

% =========================================================================
% CONFIGURATION  ← adjust this path for your machine
% =========================================================================
% Directory containing BlobFilters enhancer .m files.
blobFunctionPath = 'C:\Users\dops0035\Documents\Research\Matlab Projects\BlobFilters_sandbox\src';

% Synthetic image zoom ROI  [x_left, y_top, width, height]  (imcrop convention)
roi_synth = [30 20 70 70];

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

% =========================================================================
% SHARED PARAMETERS
% =========================================================================

% --- localThresholdFast (synthetic, 128x128, objects ~8-15 px) -----------
win_s = 21;

pSauv_s.windowsize = win_s; pSauv_s.k = 0.30; pSauv_s.r = 0.5;
pBern_s.windowsize = win_s; pBern_s.contrast = 0.05;
pBrad_s.windowsize = win_s; pBrad_s.k = 0.12;
pNibl_s.windowsize = win_s; pNibl_s.k = -0.15;

pWS_dist_s  = struct('method','distance', 'threshold',0.35, 'hMinima',1.5,  'minArea',20);
pWS_grad_s  = struct('method','gradient', 'threshold',0.35, 'smoothSigma',1.5, 'hMinima',0.03, 'minArea',20);
pWS_mark_s  = struct('method','marker',   'threshold',0.35, 'smoothSigma',1.5, 'hMinima',0.05, 'minArea',20);

% =========================================================================
% =========================================================================
%  PART A — SYNTHETIC IMAGE
% =========================================================================
% =========================================================================

fprintf('=== PART A: Synthetic image ===\n');

% -------------------------------------------------------------------------
% A1.  Generate synthetic image
% -------------------------------------------------------------------------
fprintf('Generating synthetic test image...\n');
Isynth = makeSyntheticBlobs();
fprintf('  Size: %d x %d px\n', size(Isynth,2), size(Isynth,1));

% -------------------------------------------------------------------------
% A2.  Thresholding on synthetic image
% -------------------------------------------------------------------------
fprintf('\n--- Thresholding (synthetic) ---\n');

fprintf('  sauvola...   '); tic;
BW_sauv_s = localThresholdFast(Isynth, 'method','sauvola', ...
    'windowsize',pSauv_s.windowsize, 'k',pSauv_s.k, 'r',pSauv_s.r);
fprintf('%.2fs  (%d px)\n', toc, nnz(BW_sauv_s));

fprintf('  bernsen...   '); tic;
BW_bern_s = localThresholdFast(Isynth, 'method','bernsen', ...
    'windowsize',pBern_s.windowsize, 'contrast',pBern_s.contrast);
fprintf('%.2fs  (%d px)\n', toc, nnz(BW_bern_s));

fprintf('  bradley...   '); tic;
BW_brad_s = localThresholdFast(Isynth, 'method','bradley', ...
    'windowsize',pBrad_s.windowsize, 'k',pBrad_s.k);
fprintf('%.2fs  (%d px)\n', toc, nnz(BW_brad_s));

fprintf('  niblack...   '); tic;
BW_nibl_s = localThresholdFast(Isynth, 'method','niblack', ...
    'windowsize',pNibl_s.windowsize, 'k',pNibl_s.k);
fprintf('%.2fs  (%d px)\n', toc, nnz(BW_nibl_s));

% -------------------------------------------------------------------------
% A3.  Watershed on synthetic image
% -------------------------------------------------------------------------
fprintf('\n--- Watershed (synthetic) ---\n');

fprintf('  distance...  '); tic;
[~, L_dist_s] = watershedSegment(Isynth, 'method',pWS_dist_s.method, ...
    'threshold',pWS_dist_s.threshold, 'hMinima',pWS_dist_s.hMinima, ...
    'minArea',pWS_dist_s.minArea);
fprintf('%.2fs  (%d objects)\n', toc, max(L_dist_s(:)));

fprintf('  gradient...  '); tic;
[~, L_grad_s] = watershedSegment(Isynth, 'method',pWS_grad_s.method, ...
    'threshold',pWS_grad_s.threshold, 'smoothSigma',pWS_grad_s.smoothSigma, ...
    'hMinima',pWS_grad_s.hMinima, 'minArea',pWS_grad_s.minArea);
fprintf('%.2fs  (%d objects)\n', toc, max(L_grad_s(:)));

fprintf('  marker...    '); tic;
[~, L_mark_s] = watershedSegment(Isynth, 'method',pWS_mark_s.method, ...
    'threshold',pWS_mark_s.threshold, 'smoothSigma',pWS_mark_s.smoothSigma, ...
    'hMinima',pWS_mark_s.hMinima, 'minArea',pWS_mark_s.minArea);
fprintf('%.2fs  (%d objects)\n', toc, max(L_mark_s(:)));

% -------------------------------------------------------------------------
% A4.  Figures 1-3 — synthetic
% -------------------------------------------------------------------------
rgb_sauv_s = maskRGB(BW_sauv_s);
rgb_dist_s = labelRGB(L_dist_s);
rgb_grad_s = labelRGB(L_grad_s);
rgb_mark_s = labelRGB(L_mark_s);

figure(1);
set(gcf,'Name','Synthetic — Thresholding','NumberTitle','off', ...
        'Color','k','Position',[30 680 1400 300]);
segFillFigure(1, {Isynth, single(BW_sauv_s), single(BW_bern_s), ...
                  single(BW_brad_s), single(BW_nibl_s)}, ...
    'Thresholding methods — synthetic image', ...
    {'Raw','Sauvola','Bernsen','Bradley','Niblack'}, ...
    {'gray','gray','gray','gray','gray'});

figure(2);
set(gcf,'Name','Synthetic — Watershed','NumberTitle','off', ...
        'Color','k','Position',[30 380 1400 300]);
segFillFigure(2, {Isynth, rgb_sauv_s, rgb_dist_s, rgb_grad_s, rgb_mark_s}, ...
    'Watershed methods — synthetic image', ...
    {'Raw', 'Sauvola mask', ...
     sprintf('Distance-WS (n=%d)', max(L_dist_s(:))), ...
     sprintf('Gradient-WS (n=%d)', max(L_grad_s(:))), ...
     sprintf('Marker-WS   (n=%d)', max(L_mark_s(:)))}, ...
    {'gray',[],[],[],[]});

figure(3);
set(gcf,'Name','Synthetic — Zoom marker-WS','NumberTitle','off', ...
        'Color','k','Position',[30 60 900 320]);
segFillFigure(3, ...
    {imcrop(Isynth,roi_synth), imcrop(rgb_sauv_s,roi_synth), imcrop(rgb_mark_s,roi_synth)}, ...
    sprintf('Zoom [x=%d y=%d %dx%d] — marker-WS (n=%d)', ...
            roi_synth(1),roi_synth(2),roi_synth(3),roi_synth(4),max(L_mark_s(:))), ...
    {'Raw','Sauvola mask',sprintf('Marker-WS (n=%d)',max(L_mark_s(:)))}, ...
    {'gray',[],[]});

fprintf('Figs 1-3 generated.\n');

% =========================================================================
% =========================================================================
%  PART B — REAL MITOCHONDRIAL IMAGE (CROP)
% =========================================================================
% =========================================================================

fprintf('\n=== PART B: Real mitochondrial image ===\n');

if ~hasBlobFilters
    fprintf('  BlobFilters not found at:\n    %s\n', blobFunctionPath);
    fprintf('  Cannot run Part B. Check blobFunctionPath.\n');
    fprintf('\nDone. 3 figures generated (synthetic only).\n');
    return
end

% -------------------------------------------------------------------------
% B1.  Load image
% -------------------------------------------------------------------------
mitoMat = fullfile(blobFunctionPath, '..', 'demos', 'mitImagecrop.mat');
if ~exist(mitoMat, 'file')
    fprintf('  mitImagecrop.mat not found at:\n    %s\n', mitoMat);
    fprintf('\nDone. 3 figures generated (synthetic only).\n');
    return
end
fprintf('Loading mitImagecrop.mat...\n');
tmp   = load(mitoMat);
flds  = fieldnames(tmp);
Ireal = im2single(tmp.(flds{1}));
fprintf('  Size: %d x %d px\n', size(Ireal,2), size(Ireal,1));

% -------------------------------------------------------------------------
% B2.  Cellpose ground truth
% -------------------------------------------------------------------------
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

hasCellpose = exist('cellpose','file') ~= 0;
if hasCellpose
    fprintf('  cellposeEnhance (%s, cp=%d, ft=%.1f, minSize=%d)... ', ...
            cpGTmodel, cpGTcellProb, cpGTflowThreshold, cpMinSize);
    try
        tic;
        [~, cpL] = cellposeEnhance(Ireal, pCP);
        fprintf('%.1fs  (%d objects)\n', toc, max(cpL(:)));
    catch ME
        fprintf('FAILED (%s)\n', ME.message);
        error('demoSegmentation:noCellpose', ...
              'Cellpose GT failed. Cannot continue without ground truth.');
    end
else
    error('demoSegmentation:noCellpose', ...
          'Cellpose add-on not installed. Cannot generate ground truth.');
end

cpBW  = logical(cpL > 0);
cpRGB = labelRGB(cpL);
nCP   = max(cpL(:));

% -------------------------------------------------------------------------
% B3.  OGS preprocessing
% -------------------------------------------------------------------------
fprintf('\n=== OGS preprocessing ===\n');
pOGS.sigmaAlong = 4;  pOGS.sigmaAcross = 1.5;
pOGS.orientations = 8;  pOGS.sigmaGrad = 1.5;  pOGS.sigmaInt = 5;
fprintf('  orientedGaussSmooth... '); tic;
Ism = orientedGaussSmooth(Ireal, pOGS);
fprintf('%.1fs\n', toc);

% -------------------------------------------------------------------------
% B4.  Watershed sweep grids (per-task threshold ranges)
% -------------------------------------------------------------------------
hMinimaVals = [0.02 0.05 0.10];
minAreaVals = [20 50];

% Default: log, fiber, rod
wsThreshDef     = [0.20 0.25 0.30 0.35 0.40];
wsThreshLowDef  = {[], 0.10, 0.15};

% Capsule: higher thresholds, raised hysteresis lower bounds
wsThreshHigh    = [0.30 0.35 0.40 0.45 0.50];
wsThreshLowHigh = {[], 0.20, 0.25};

% -------------------------------------------------------------------------
% B5.  Enhancer task definitions + per-task WS config
% -------------------------------------------------------------------------
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

% -------------------------------------------------------------------------
% B6.  Per-enhancer sweep — one figure per enhancer (Figs 4..4+nTasks-1)
% -------------------------------------------------------------------------
fprintf('\n=== Enhancer parameter sweep ===\n');
gtLabel = sprintf('%s cp=%d ft=%.1f minSize=%d', ...
                  cpGTmodel, cpGTcellProb, cpGTflowThreshold, cpMinSize);

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
    bestDescPerCfg   = cell(1, nCfgs);   % panel subtitle (cfgDesc + F1)

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

    % Figure: Raw | cfg1-bestWS | cfg2-bestWS | ... | Cellpose GT
    panels      = [{Ireal}, cellfun(@labelRGB, bestLabelsPerCfg, 'UniformOutput',false), {cpRGB}];
    panelTitles = [{'Raw'}, bestDescPerCfg, {sprintf('Cellpose GT\n(n=%d)', nCP)}];
    cmaps       = [{'gray'}, repmat({[]}, 1, nCfgs), {[]}];

    figNum = 3 + ti;
    figW   = max(900, 300 * (nCfgs + 2));
    figure(figNum);
    set(gcf, 'Name', sprintf('%s sweep', tk.name), 'NumberTitle','off', ...
             'Color','k', 'Position',[30, 50, figW, 360]);
    segFillFigure(figNum, panels, ...
        sprintf('%s — best WS per config  |  best: %s, F1=%.3f  |  GT: %s', ...
                tk.name, tk.cfgDescs{iBest}, bestF1overall, gtLabel), ...
        panelTitles, cmaps);

    fprintf('  >> Best overall: config %d (%s)  F1=%.4f\n', ...
            iBest, tk.cfgDescs{iBest}, bestF1overall);
    fprintf('  Fig %d generated.\n', figNum);
end

fprintf('\nDone. %d figures generated.\n', 3 + nTasks);


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function I = makeSyntheticBlobs()
% makeSyntheticBlobs  128x128 synthetic image for segmentation demos.
%   Four structures on a noisy background:
%     1+2. Two touching compact blobs — watershed should separate these.
%     3.   Elongated horizontal rod.
%     4.   Small punctum (below default minArea=20; tests area filtering).
[xx, yy] = meshgrid(1:128, 1:128);
rng(42);
blob1 = 0.90 * exp(-((xx-50).^2/18 + (yy-50).^2/18));
blob2 = 0.85 * exp(-((xx-70).^2/18 + (yy-50).^2/18));
rod   = 0.75 * exp(-((xx-64).^2/200 + (yy-90).^2/10));
dot   = 0.60 * exp(-((xx-105).^2/4  + (yy-105).^2/4));
I = blob1 + blob2 + rod + dot;
I = I / max(I(:));
I = I + 0.05 * randn(128, 128, 'single');
I = im2single(max(0, min(1, I)));
end


function rgb = labelRGB(L)
% labelRGB  Label image → uint8 RGB with black background.
if max(L(:)) == 0
    rgb = zeros(size(L,1), size(L,2), 3, 'uint8');
else
    rgb = label2rgb(L, 'jet', 'k', 'shuffle');
end
end


function rgb = maskRGB(BW)
% maskRGB  Logical mask → uint8 RGB greyscale (for truecolor panels).
bw8 = uint8(BW) * 255;
rgb = cat(3, bw8, bw8, bw8);
end


function segFillFigure(figNum, panels, figTitle, panelTitles, cmaps)
% segFillFigure  Populate a 1xN figure with black background.
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


function tasks = buildTasks()
% buildTasks  Define all enhancer tasks with their parameter configuration grids.

tasks = {};

% --- logEnhance -----------------------------------------------------------
t.name = 'logEnhance';  t.tag = 'log-WS';
t.fn   = @(I, p) logEnhance(I, p);
t.cfgs = {
    struct('sigmas',[1 2 3],      'normalize',true),
    struct('sigmas',[2 3 4 5],    'normalize',true),
    struct('sigmas',[2 3 4 5 6],  'normalize',true),
    struct('sigmas',[3 4 5 6 7],  'normalize',true)
};
t.cfgDescs = {
    'sigmas=[1 2 3]',
    'sigmas=[2 3 4 5]',
    'sigmas=[2 3 4 5 6]',
    'sigmas=[3 4 5 6 7]'
};
t.maxExpand = 4;
tasks{end+1} = t;

% --- fiberEnhance ---------------------------------------------------------
t.name = 'fiberEnhance';  t.tag = 'fib-WS';
t.fn   = @(I, p) fiberEnhance(I, p);
t.cfgs = {
    struct('widths',[4 5 6 7],    'multimode','stack','normalize',true),
    struct('widths',[6 7 8 9 10], 'multimode','stack','normalize',true),
    struct('widths',[8 9 10 12],  'multimode','stack','normalize',true)
};
t.cfgDescs = {
    'widths=[4 5 6 7]',
    'widths=[6 7 8 9 10]',
    'widths=[8 9 10 12]'
};
t.maxExpand = 7;
tasks{end+1} = t;

% --- capsuleEnhance single ------------------------------------------------
t.name = 'capsule(single)';  t.tag = 'capS-WS';
t.fn   = @(I, p) capsuleEnhance(I, p);
t.cfgs = {
    struct('lengths',[8 12 16],      'width',6,'orientations',8, 'mode','single','normalize',true),
    struct('lengths',[12 16 20],     'width',8,'orientations',8, 'mode','single','normalize',true),
    struct('lengths',[12 16 20 28],  'width',8,'orientations',8, 'mode','single','normalize',true)
};
t.cfgDescs = {
    'lens=[8 12 16] w=6 or=8',
    'lens=[12 16 20] w=8 or=8',
    'lens=[12 16 20 28] w=8 or=8'
};
t.maxExpand = 5;
tasks{end+1} = t;

% --- capsuleEnhance DoC ---------------------------------------------------
t.name = 'capsule(DoC)';  t.tag = 'capD-WS';
t.fn   = @(I, p) capsuleEnhance(I, p);
t.cfgs = {
    struct('lengths',[8 12 16 20],       'width',6,'wideWidth',14,'alpha',0.55,'orientations',8,'mode','doc','normalize',true),
    struct('lengths',[12 16 20 28],      'width',8,'wideWidth',18,'alpha',0.55,'orientations',8,'mode','doc','normalize',true),
    struct('lengths',[12 16 20 28 36 40],'width',8,'wideWidth',18,'alpha',0.55,'orientations',8,'mode','doc','normalize',true)
};
t.cfgDescs = {
    'lens=[8..20] w=6 wW=14 or=8',
    'lens=[12..28] w=8 wW=18 or=8',
    'lens=[12..40] w=8 wW=18 or=8'
};
t.maxExpand = 5;
tasks{end+1} = t;

% --- rodGranulometryEnhance -----------------------------------------------
t.name = 'rodGranulometry';  t.tag = 'rod-WS';
t.fn   = @(I, p) rodGranulometryEnhance(I, p);
t.cfgs = {
    struct('lengths',[4 8 12],          'orientations',8, 'normalize',true),
    struct('lengths',[8 12 16 20],      'orientations',8, 'normalize',true),
    struct('lengths',[8 12 16 20 28 36],'orientations',8, 'normalize',true)
};
t.cfgDescs = {
    'lens=[4 8 12] or=8',
    'lens=[8 12 16 20] or=8',
    'lens=[8 12 16 20 28 36] or=8'
};
t.maxExpand = 8;
tasks{end+1} = t;
end


function [f1, prec, rec, iou] = evalBinary(BW, GT)
% evalBinary  Pixel-level binary segmentation metrics vs ground truth.
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
% wsDesc  Short string describing a watershed parameter set.
if isempty(ws.threshLow)
    s = sprintf('thr=%.2f h=%.2f minA=%d', ws.threshold, ws.hMinima, ws.minArea);
else
    s = sprintf('thr=%.2f lo=%.2f h=%.2f minA=%d', ...
                ws.threshold, ws.threshLow, ws.hMinima, ws.minArea);
end
end
