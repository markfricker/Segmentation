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
%   Fig 4 — capsule(DoC) sweep
%   Fig 5 — rodGranulometry sweep
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

    % Figure: Raw | cfg1-bestWS | cfg2-bestWS | ... | Cellpose GT
    panels      = [{Ireal}, cellfun(@labelRGB, bestLabelsPerCfg, 'UniformOutput',false), {cpRGB}];
    panelTitles = [{'Raw'}, bestDescPerCfg, {sprintf('Cellpose GT\n(n=%d)', nCP)}];
    cmaps       = [{'gray'}, repmat({[]}, 1, nCfgs), {[]}];

    figNum = ti;
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

fprintf('\nDone. %d figures generated.\n', nTasks);


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


function tasks = buildTasks()
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
    struct('lengths',[8 12 16],      'width',6,'orientations',8,'mode','single','normalize',true),
    struct('lengths',[12 16 20],     'width',8,'orientations',8,'mode','single','normalize',true),
    struct('lengths',[12 16 20 28],  'width',8,'orientations',8,'mode','single','normalize',true)
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
    struct('lengths',[4 8 12],          'orientations',8,'normalize',true),
    struct('lengths',[8 12 16 20],      'orientations',8,'normalize',true),
    struct('lengths',[8 12 16 20 28 36],'orientations',8,'normalize',true)
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
