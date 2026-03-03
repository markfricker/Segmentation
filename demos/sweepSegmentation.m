%%
%%
% sweepSegmentation.m
%
% Parameter sweep: find optimal enhancer + watershed configuration for
% mitochondria segmentation, benchmarked against Cellpose as ground truth.
% Includes a boundary-refinement stage (refineSegment) applied to each
% best watershed result to compare Chan-Vese AC vs Voronoi-dilate methods.
%
% STRATEGY
%   For each enhancer (logEnhance, fiberEnhance, capsule-single,
%   capsule-DoC, rodGranulometryEnhance):
%     - Test N enhancer parameter configurations (scale ranges, widths, etc.)
%     - For each configuration, sweep marker-watershed parameters
%       (threshold, hMinima, minArea)
%     - Evaluate binary segmentation quality vs Cellpose GT at every combo
%     - Keep the best-F1 result for each enhancer
%   Then, for each best watershed result:
%     - Apply refineSegment with 'chanvese' and 'dilate' methods
%     - Report F1/Prec/Rec/IoU delta vs the unrefined watershed result
%
% METRICS  (pixel-level, binary: watershed foreground vs Cellpose foreground)
%   F1 / Dice  =  2|WS∩CP| / (|WS|+|CP|)              <- primary ranking
%   Precision  =  |WS∩CP| / |WS|
%   Recall     =  |WS∩CP| / |CP|
%   IoU        =  |WS∩CP| / |WS∪CP|
%
% FIGURES PRODUCED
%   Fig 1 — Best watershed result per enhancer vs Cellpose (full FOV)
%           Raw | log-WS | fib-WS | capS-WS | capD-WS | rod-WS | CP
%   Fig 2 — Zoom cluster 1   (same 7 columns)
%   Fig 3 — Zoom cluster 2   (same 7 columns)
%   Fig 4 — Refinement comparison — full FOV
%           Raw | best-WS | +chanvese | +dilate | Cellpose GT
%   Fig 5 — Refinement comparison — zoom cluster 1
%   Fig 6 — Refinement comparison — zoom cluster 2
%
% CONSOLE OUTPUT
%   Ranked table of best watershed parameters per enhancer.
%   Refinement comparison table: WS → chanvese → dilate metrics.
%
% REQUIREMENTS
%   localThresholdFast.m, watershedSegment.m, refineSegment.m on path.
%   BlobFilters enhancers and cellposeEnhance on path (set blobFunctionPath).
%   Cellpose add-on + Python cellpose for ground-truth generation.
%   Image Processing Toolbox (including activecontour for Chan-Vese).
%
% NOTE ON RUNTIME
%   Enhancement computation dominates (~2-4 min total).
%   Watershed sweep (480 calls) adds ~1 min.
%   Chan-Vese refinement: ~5-30 s per enhancer depending on object count.
%   Total: ~5-8 min.  If cpL already exists in the workspace (from
%   demoSegmentation), Cellpose will not be re-run.

clear; clc; close all;

% =========================================================================
% CONFIGURATION
% =========================================================================
blobFunctionPath = 'C:\Users\dops0035\Documents\Research\Matlab Projects\BlobFilters_sandbox\src';
blobImagePath    = 'C:\Users\dops0035\Documents\Research\Matlab Projects\BlobFilters_sandbox\src\.claude\worktrees\thirsty-wescoff\demos';

roi1 = [100  80  220 220];   % zoom cluster 1  [x y w h]
roi2 = [330 220  220 220];   % zoom cluster 2

% ---- refinement parameters (section 9) ----------------------------------
% maxExpand is set per-enhancer inside buildTasks() — tune it there.
refNIter       = 100;   % Chan-Vese AC iterations
refFgThresh    = 0.10;  % dilate method: minimum raw intensity to accept
refMinArea     = 20;    % minimum object area after refinement (pixels)

% =========================================================================
% Path setup
% =========================================================================
demoDir = fileparts(mfilename('fullpath'));
addpath(fullfile(demoDir, '..', 'src'));

if ~isempty(blobFunctionPath) && exist(blobFunctionPath, 'dir')
    addpath(blobFunctionPath);
end
hasBlobFilters = exist('rodGranulometryEnhance', 'file') ~= 0;
if ~hasBlobFilters
    error('sweepSegmentation:noBlobFilters', ...
          'BlobFilters enhancers not found. Set blobFunctionPath correctly.');
end

% =========================================================================
% 1.  Load real mitochondrial image
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
    error('sweepSegmentation:noImage', ...
          'Image not found. Check blobImagePath.');
end

% =========================================================================
% 2.  Cellpose ground truth
%     Reuse from workspace if already computed (e.g. from demoSegmentation).
% =========================================================================
fprintf('\n=== Cellpose ground truth ===\n');
hasCellpose = exist('cellpose','file') ~= 0;

if evalin('base','exist(''cpL'',''var'')') && ...
        isequal(size(evalin('base','cpL')), size(Ireal))
    cpL  = evalin('base', 'cpL');
    fprintf('  Reusing cpL from workspace  (n=%d objects)\n', max(cpL(:)));
elseif hasCellpose
    pCP.model = 'cyto3';  pCP.diameter = 10;
    pCP.cellProb = -2;    pCP.flowThreshold = 0.8;
    fprintf('  Running cellposeEnhance...  ');
    try
        tic;
        [~, cpL] = cellposeEnhance(Ireal, pCP);
        fprintf('%.1fs  (n=%d objects)\n', toc, max(cpL(:)));
    catch ME
        error('sweepSegmentation:cellposeFailed', ...
              'Cellpose failed: %s', ME.message);
    end
else
    error('sweepSegmentation:noCellpose', ...
          'Cellpose add-on not installed. Cannot generate ground truth.');
end

cpBW = logical(cpL > 0);   % binary ground-truth mask
cpRGB = labelRGB(cpL);
nCP   = max(cpL(:));

% =========================================================================
% 3.  OGS preprocessing
% =========================================================================
fprintf('\n=== OGS preprocessing ===\n');
pOGS.sigmaAlong = 4;  pOGS.sigmaAcross = 1.5;
pOGS.orientations = 8;  pOGS.sigmaGrad = 1.5;  pOGS.sigmaInt = 5;
fprintf('  orientedGaussSmooth...  '); tic;
Ism = orientedGaussSmooth(Ireal, pOGS);
fprintf('%.1fs\n', toc);

% =========================================================================
% 4.  Watershed sweep grid
% =========================================================================
thresholds    = [0.20 0.25 0.30 0.35 0.40];
hMinimaVals   = [0.02 0.05 0.10];
minAreaVals   = [20 50];
threshLowVals = {[], 0.10, 0.15};   % [] = single-threshold; values = hysteresis low

wsGrid = {};
for thr = thresholds
    for hm = hMinimaVals
        for ma = minAreaVals
            for tli = 1:numel(threshLowVals)
                tl = threshLowVals{tli};
                if ~isempty(tl) && tl >= thr, continue; end   % low must be < high
                wsGrid{end+1} = struct('threshold',thr,'hMinima',hm, ...
                                       'minArea',ma,'threshLow',tl); %#ok<SAGROW>
            end
        end
    end
end
nWS = numel(wsGrid);

% =========================================================================
% 5.  Enhancer task definitions
%
%     Each task: .name  (display string)
%                .tag   (short identifier for output)
%                .fn    (function handle: fn(Ism, cfg) -> R)
%                .cfgs  (cell array of parameter structs)
%                .cfgDescs (cell array of description strings)
% =========================================================================

tasks = buildTasks();
nTasks = numel(tasks);

% =========================================================================
% 6.  Main sweep
% =========================================================================
fprintf('\n=== Parameter sweep (%d enhancer configs x %d WS combos) ===\n', ...
        sum(cellfun(@(t) numel(t.cfgs), tasks)), nWS);

results    = cell(nTasks, 1);   % best result struct per task
bestLabels = cell(nTasks, 1);   % best label image per task
bestEnhMap = cell(nTasks, 1);   % best enhancement map per task
resultsHyst = cell(nTasks, 1);  % best hysteresis-only result per task
labelsHyst  = cell(nTasks, 1);  % best hysteresis label image per task
resultsStd  = cell(nTasks, 1);  % best standard-threshold result per task
labelsStd   = cell(nTasks, 1);  % best standard-threshold label image per task

for ti = 1:nTasks
    tk = tasks{ti};
    fprintf('\n--- %s (%d configs x %d WS = %d evals) ---\n', ...
            tk.name, numel(tk.cfgs), nWS, numel(tk.cfgs)*nWS);

    emptyRes = struct('f1',0,'prec',0,'rec',0,'iou',0, ...
                      'enhIdx',1,'wsIdx',1,'enhDesc','','wsDesc','');
    best     = emptyRes;
    bestHyst = emptyRes;
    bestStd  = emptyRes;
    bestL     = zeros(size(Ireal), 'uint16');
    bestR     = zeros(size(Ireal), 'single');
    bestHystL = zeros(size(Ireal), 'uint16');
    bestStdL  = zeros(size(Ireal), 'uint16');

    for ci = 1:numel(tk.cfgs)
        fprintf('  Config %d/%d: %s\n    Computing enhancement... ', ...
                ci, numel(tk.cfgs), tk.cfgDescs{ci});
        tic;
        try
            R = tk.fn(Ism, tk.cfgs{ci});
            fprintf('%.1fs\n', toc);
        catch ME
            fprintf('SKIPPED (%s)\n', ME.message);
            continue;
        end

        fprintf('    Sweeping %d watershed combos: ', nWS);
        t0 = tic;
        for wi = 1:nWS
            ws = wsGrid{wi};
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

            [f1, prec, rec, iou] = evalBinary(BW, cpBW);

            newRes = struct('f1',f1,'prec',prec,'rec',rec,'iou',iou, ...
                            'enhIdx',ci,'wsIdx',wi, ...
                            'enhDesc',tk.cfgDescs{ci},'wsDesc',wsDesc(ws));
            if f1 > best.f1
                best  = newRes;
                bestL = L;
                bestR = R;
            end
            if isempty(ws.threshLow) && f1 > bestStd.f1
                bestStd  = newRes;
                bestStdL = L;
            end
            if ~isempty(ws.threshLow) && f1 > bestHyst.f1
                bestHyst  = newRes;
                bestHystL = L;
            end
        end
        fprintf('done (%.1fs)\n', toc(t0));
    end

    results{ti}     = best;
    bestLabels{ti}  = bestL;
    bestEnhMap{ti}  = bestR;
    resultsHyst{ti} = bestHyst;
    labelsHyst{ti}  = bestHystL;
    resultsStd{ti}  = bestStd;
    labelsStd{ti}   = bestStdL;

    fprintf('  >> Best F1=%.4f  Prec=%.4f  Rec=%.4f  IoU=%.4f\n', ...
            best.f1, best.prec, best.rec, best.iou);
    fprintf('     Enhancer: %s\n', best.enhDesc);
    fprintf('     Watershed: %s\n', best.wsDesc);
end

% =========================================================================
% 7.  Summary table
% =========================================================================
fprintf('\n');
printTable(tasks, results, nCP);

% =========================================================================
% 8.  Figures
% =========================================================================
panelTitles = {'Raw'};
for ti = 1:nTasks
    panelTitles{end+1} = sprintf('%s\nF1=%.3f n=%d', ...
        tasks{ti}.tag, results{ti}.f1, max(bestLabels{ti}(:))); %#ok<SAGROW>
end
panelTitles{end+1} = sprintf('Cellpose GT\n(n=%d)', nCP);
panelCmaps = [{'gray'}, repmat({[]}, 1, nTasks), {[]}];

% bestLabels is Nx1; transpose to row before cellfun so horzcat works
panels = [{Ireal}, cellfun(@labelRGB, bestLabels(:)', 'UniformOutput', false), {cpRGB}];

figure(1);
set(gcf,'Name','Best result per enhancer — full FOV','NumberTitle','off', ...
        'Color','k','Position',[30 680 1820 340]);
segFillFigure(1, panels, ...
    sprintf('Best marker-WS per enhancer vs Cellpose GT (n=%d) — full FOV', nCP), ...
    panelTitles, panelCmaps);

figure(2);
set(gcf,'Name','Best result per enhancer — zoom 1','NumberTitle','off', ...
        'Color','k','Position',[30 340 1820 340]);
segFillFigure(2, cellfun(@(p) imcrop(p,roi1), panels, 'UniformOutput',false), ...
    sprintf('Zoom cluster 1  [x=%d y=%d %dx%d]', roi1(1),roi1(2),roi1(3),roi1(4)), ...
    panelTitles, panelCmaps);

figure(3);
set(gcf,'Name','Best result per enhancer — zoom 2','NumberTitle','off', ...
        'Color','k','Position',[30 0 1820 340]);
segFillFigure(3, cellfun(@(p) imcrop(p,roi2), panels, 'UniformOutput',false), ...
    sprintf('Zoom cluster 2  [x=%d y=%d %dx%d]', roi2(1),roi2(2),roi2(3),roi2(4)), ...
    panelTitles, panelCmaps);

fprintf('\nDone. 3 figures generated.\n');

% =========================================================================
% 8b.  Hysteresis vs standard comparison
% =========================================================================
fprintf('\n');
printHystTable(tasks, resultsStd, resultsHyst, nCP);

% Figures 7-9: best hysteresis result per enhancer (parallel to Figs 1-3)
panelTitlesH = {'Raw'};
for ti = 1:nTasks
    if resultsHyst{ti}.f1 > 0
        lbl = sprintf('%s\nF1=%.3f n=%d', tasks{ti}.tag, ...
                      resultsHyst{ti}.f1, max(labelsHyst{ti}(:)));
    else
        lbl = sprintf('%s\n(no hyst result)', tasks{ti}.tag);
    end
    panelTitlesH{end+1} = lbl; %#ok<SAGROW>
end
panelTitlesH{end+1} = sprintf('Cellpose GT\n(n=%d)', nCP);
panelCmapsH = [{'gray'}, repmat({[]}, 1, nTasks), {[]}];
panelsH = [{Ireal}, cellfun(@labelRGB, labelsHyst(:)', 'UniformOutput', false), {cpRGB}];

figure(7);
set(gcf,'Name','Best hysteresis result per enhancer — full FOV','NumberTitle','off', ...
        'Color','k','Position',[30 680 1820 340]);
segFillFigure(7, panelsH, ...
    sprintf('Best hysteresis-WS per enhancer vs Cellpose GT (n=%d) — full FOV', nCP), ...
    panelTitlesH, panelCmapsH);

figure(8);
set(gcf,'Name','Best hysteresis result per enhancer — zoom 1','NumberTitle','off', ...
        'Color','k','Position',[30 340 1820 340]);
segFillFigure(8, cellfun(@(p) imcrop(p,roi1), panelsH, 'UniformOutput',false), ...
    sprintf('Zoom cluster 1  [x=%d y=%d %dx%d]', roi1(1),roi1(2),roi1(3),roi1(4)), ...
    panelTitlesH, panelCmapsH);

figure(9);
set(gcf,'Name','Best hysteresis result per enhancer — zoom 2','NumberTitle','off', ...
        'Color','k','Position',[30 0 1820 340]);
segFillFigure(9, cellfun(@(p) imcrop(p,roi2), panelsH, 'UniformOutput',false), ...
    sprintf('Zoom cluster 2  [x=%d y=%d %dx%d]', roi2(1),roi2(2),roi2(3),roi2(4)), ...
    panelTitlesH, panelCmapsH);

fprintf('Hysteresis figures done. Figs 7-9 generated.\n');

% =========================================================================
% 9.  Boundary refinement: apply refineSegment to each best WS result
% =========================================================================
fprintf('\n=== Boundary refinement ===\n');
fprintf('    maxExpand: per-enhancer (see buildTasks)  |  chanvese: nIter=%d  |  dilate: fgThr=%.2f\n', ...
        refNIter, refFgThresh);

refMethods   = {'chanvese', 'dilate'};
nRefMethods  = numel(refMethods);
resultsRef   = cell(nTasks, nRefMethods);   % metrics structs
refinedL     = cell(nTasks, nRefMethods);   % refined label images

for ti = 1:nTasks
    L_ws = bestLabels{ti};
    for mi = 1:nRefMethods
        mName = refMethods{mi};
        fprintf('  %-16s + %-8s ... ', tasks{ti}.tag, mName);
        tic;
        try
            [BW_r, L_r] = refineSegment(Ireal, L_ws, ...
                'method',      mName, ...
                'maxExpand',   tasks{ti}.maxExpand, ...
                'nIter',       refNIter, ...
                'fgThreshold', refFgThresh, ...
                'minArea',     refMinArea);
            [f1r, pr, rr, iour] = evalBinary(BW_r, cpBW);
            resultsRef{ti, mi} = struct('f1',f1r,'prec',pr,'rec',rr,'iou',iour);
            refinedL{ti, mi}   = L_r;
            dF1 = f1r - results{ti}.f1;
            fprintf('%.1fs  F1=%.4f  (%+.4f vs WS)\n', toc, f1r, dF1);
        catch ME
            fprintf('SKIPPED (%s)\n', ME.message);
            resultsRef{ti, mi} = struct('f1',0,'prec',0,'rec',0,'iou',0);
            refinedL{ti, mi}   = zeros(size(Ireal), 'uint16');
        end
    end
end

% =========================================================================
% 10. Refinement comparison table
% =========================================================================
fprintf('\n');
printRefinementTable(tasks, results, resultsRef, refMethods, nCP);

% =========================================================================
% 11. Refinement figures — best overall WS result + both refinements
%     (select the enhancer with the highest post-refinement chanvese F1)
% =========================================================================
% Find which task benefits most from chanvese refinement
cvF1s    = cellfun(@(r) r.f1, resultsRef(:, 1));   % chanvese column
[~, iBest] = max(cvF1s);

bestTag   = tasks{iBest}.tag;
L_ws_best = bestLabels{iBest};
L_cv_best = refinedL{iBest, 1};
L_dil_best = refinedL{iBest, 2};
wsF1  = results{iBest}.f1;
cvF1  = resultsRef{iBest, 1}.f1;
dilF1 = resultsRef{iBest, 2}.f1;

refPanelTitles = { ...
    'Raw', ...
    sprintf('%s  WS\nF1=%.3f', bestTag, wsF1), ...
    sprintf('+chanvese\nF1=%.3f  (%+.3f)', cvF1,  cvF1  - wsF1), ...
    sprintf('+dilate\nF1=%.3f  (%+.3f)',   dilF1, dilF1 - wsF1), ...
    sprintf('Cellpose GT\n(n=%d)', nCP) };
refCmaps = {'gray', [], [], [], []};

refPanels = { ...
    Ireal, ...
    labelRGB(L_ws_best), ...
    labelRGB(L_cv_best), ...
    labelRGB(L_dil_best), ...
    cpRGB };

figure(4);
set(gcf,'Name','Refinement comparison — full FOV','NumberTitle','off', ...
        'Color','k','Position',[30 680 1820 340]);
segFillFigure(4, refPanels, ...
    sprintf('Refinement  (%s): WS → Chan-Vese / Voronoi-dilate vs Cellpose GT', bestTag), ...
    refPanelTitles, refCmaps);

figure(5);
set(gcf,'Name','Refinement comparison — zoom 1','NumberTitle','off', ...
        'Color','k','Position',[30 340 1820 340]);
segFillFigure(5, cellfun(@(p) imcrop(p, roi1), refPanels, 'UniformOutput', false), ...
    sprintf('Zoom cluster 1  [x=%d y=%d %dx%d]', roi1(1),roi1(2),roi1(3),roi1(4)), ...
    refPanelTitles, refCmaps);

figure(6);
set(gcf,'Name','Refinement comparison — zoom 2','NumberTitle','off', ...
        'Color','k','Position',[30 0 1820 340]);
segFillFigure(6, cellfun(@(p) imcrop(p, roi2), refPanels, 'UniformOutput', false), ...
    sprintf('Zoom cluster 2  [x=%d y=%d %dx%d]', roi2(1),roi2(2),roi2(3),roi2(4)), ...
    refPanelTitles, refCmaps);

fprintf('\nDone. 6 figures generated.\n');


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function tasks = buildTasks()
% Define all enhancer tasks with their parameter configuration grids.

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
t.maxExpand = 4;   % LoG blobs are compact; conservative expansion
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
t.maxExpand = 7;   % fiber boundaries fade gradually; allow more expansion
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
t.maxExpand = 5;   % capsule width sets the scale; moderate expansion
tasks{end+1} = t;

% --- capsuleEnhance DoC ---------------------------------------------------
t.name = 'capsule(DoC)';  t.tag = 'capD-WS';
t.fn   = @(I, p) capsuleEnhance(I, p);
t.cfgs = {
    struct('lengths',[8 12 16 20],      'width',6,'wideWidth',14,'alpha',0.55,'orientations',8, 'mode','doc','normalize',true),
    struct('lengths',[12 16 20 28],     'width',8,'wideWidth',18,'alpha',0.55,'orientations',8, 'mode','doc','normalize',true),
    struct('lengths',[12 16 20 28 36 40],'width',8,'wideWidth',18,'alpha',0.55,'orientations',8, 'mode','doc','normalize',true)
};
t.cfgDescs = {
    'lens=[8..20] w=6 wW=14 or=8',
    'lens=[12..28] w=8 wW=18 or=8',
    'lens=[12..40] w=8 wW=18 or=8'
};
t.maxExpand = 5;   % DoC suppresses background well; same scale as single
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
t.maxExpand = 8;   % rod halos are diffuse; allow largest expansion
tasks{end+1} = t;
end


function printHystTable(tasks, resultsStd, resultsHyst, nCP)
% printHystTable  Print standard vs hysteresis best-F1 comparison per enhancer.
nT = numel(tasks);

w1  = 16;
sep = repmat('=', 1, w1 + 60);
fprintf('%s\n', sep);
fprintf('HYSTERESIS vs STANDARD THRESHOLD — best F1 per enhancer (Cellpose GT n=%d)\n', nCP);
fprintf('%s\n', sep);
fprintf(' %-*s  %8s  %8s  %8s  %s\n', w1, 'Enhancer', 'Std F1', 'Hyst F1', 'Delta', 'Best hyst params');
fprintf(' %s\n', repmat('-', 1, numel(sep)-1));
for ti = 1:nT
    rs = resultsStd{ti};
    rh = resultsHyst{ti};
    if rh.f1 > 0
        dF1 = rh.f1 - rs.f1;
        sign_str = '+';
        if dF1 < 0, sign_str = ''; end
        fprintf(' %-*s  %8.4f  %8.4f  %s%.4f  %s\n', ...
                w1, tasks{ti}.name, rs.f1, rh.f1, sign_str, dF1, rh.wsDesc);
    else
        fprintf(' %-*s  %8.4f  %8s  %8s  (no hysteresis combos ran)\n', ...
                w1, tasks{ti}.name, rs.f1, '—', '—');
    end
end
fprintf(' %s\n', repmat('-', 1, numel(sep)-1));
fprintf('%s\n\n', sep);
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


function printTable(tasks, results, nCP)
% printTable  Print a ranked summary table to the console.
nT = numel(tasks);

% Sort by F1 descending
f1s  = cellfun(@(r) r.f1,  results);
[~, ord] = sort(f1s, 'descend');

w1 = 18;  w2 = 30;  w3 = 22;  % column widths
sep = repmat('=', 1, w1+w2+w3+30);

fprintf('%s\n', sep);
fprintf('PARAMETER SWEEP — marker-watershed vs Cellpose GT (n=%d)\n', nCP);
fprintf('%s\n', sep);
fprintf(' %-4s  %-*s  %-*s  %-*s  %6s  %6s  %6s  %6s\n', ...
        'Rank', w1,'Enhancer', w2,'Best enhancer config', w3,'Watershed params', ...
        'F1', 'Prec', 'Rec', 'IoU');
fprintf(' %s\n', repmat('-',1,numel(sep)-1));

for k = 1:nT
    ti = ord(k);
    r  = results{ti};
    fprintf(' %-4d  %-*s  %-*s  %-*s  %6.4f  %6.4f  %6.4f  %6.4f\n', ...
            k, w1, tasks{ti}.name, w2, r.enhDesc, w3, r.wsDesc, ...
            r.f1, r.prec, r.rec, r.iou);
end

fprintf(' %s\n', repmat('-',1,numel(sep)-1));
fprintf(' %-4s  %-*s  %-*s  %-*s  %6s  %6s  %6s  %6s\n', ...
        'REF', w1,'Cellpose', w2, ...
        'model=cyto3 d=10 prob=-2 flow=0.8', w3, '(ground truth)', ...
        '1.0000','1.0000','1.0000','1.0000');
fprintf('%s\n\n', sep);

% Additional note on best classical vs Cellpose gap
bestF1 = max(f1s);
fprintf('Best classical F1 = %.4f   Gap to Cellpose = %.4f (%.1f%%)\n', ...
        bestF1, 1 - bestF1, (1 - bestF1)*100);
end


function printRefinementTable(tasks, results, resultsRef, refMethods, nCP)
% printRefinementTable  Print WS vs refined-methods comparison table.
nT  = numel(tasks);
nM  = numel(refMethods);

w1 = 16;   % enhancer name column width

% Sort by best post-refinement F1 (max over all methods including WS)
allF1 = zeros(nT, 1);
for ti = 1:nT
    f1s = results{ti}.f1;
    for mi = 1:nM
        f1s = max(f1s, resultsRef{ti, mi}.f1);
    end
    allF1(ti) = f1s;
end
[~, ord] = sort(allF1, 'descend');

sep = repmat('=', 1, w1 + 8 + (nM+1)*28 + 4);

fprintf('%s\n', sep);
fprintf('REFINEMENT COMPARISON — marker-WS + refineSegment vs Cellpose GT (n=%d)\n', nCP);
fprintf('    Metrics: F1 / Prec / Rec\n');
fprintf('%s\n', sep);

% Header
hdr = sprintf(' %-*s  %26s', w1, 'Enhancer', 'Watershed (no refinement)');
for mi = 1:nM
    hdr = [hdr, sprintf('  %26s', refMethods{mi})]; %#ok<AGROW>
end
fprintf('%s\n', hdr);
fprintf(' %s\n', repmat('-', 1, numel(sep)-1));

for k = 1:nT
    ti = ord(k);
    r  = results{ti};
    line = sprintf(' %-*s  F1=%6.4f P=%6.4f R=%6.4f', ...
                   w1, tasks{ti}.name, r.f1, r.prec, r.rec);
    for mi = 1:nM
        rr = resultsRef{ti, mi};
        dF1 = rr.f1 - r.f1;
        line = [line, sprintf('  F1=%6.4f P=%6.4f (%+.4f)', ...
                              rr.f1, rr.prec, dF1)]; %#ok<AGROW>
    end
    fprintf('%s\n', line);
end

fprintf(' %s\n', repmat('-', 1, numel(sep)-1));

% Best overall
bestWS = max(cellfun(@(r) r.f1, results));
bestCV = max(max(cellfun(@(r) r.f1, resultsRef(:,1))));
if nM > 1
    bestDil = max(max(cellfun(@(r) r.f1, resultsRef(:,2))));
else
    bestDil = NaN;
end
fprintf('\n Best WS F1 = %.4f', bestWS);
fprintf('   Best +chanvese F1 = %.4f  (%+.4f)', bestCV,  bestCV  - bestWS);
if ~isnan(bestDil)
    fprintf('   Best +dilate F1 = %.4f  (%+.4f)', bestDil, bestDil - bestWS);
end
fprintf('   Cellpose = 1.0000\n');
fprintf('%s\n\n', sep);
end


function rgb = labelRGB(L)
% labelRGB  Label image -> uint8 RGB, black background.
if max(L(:)) == 0
    rgb = zeros(size(L,1), size(L,2), 3, 'uint8');
else
    rgb = label2rgb(L, 'jet', 'k', 'shuffle');
end
end


function segFillFigure(figNum, panels, figTitle, panelTitles, panelCmaps)
% segFillFigure  Populate a 1xN figure with black background.
figure(figNum);
N = numel(panels);
for k = 1:N
    ax = subplot(1, N, k);
    if isempty(panelCmaps{k})
        imshow(panels{k});
    else
        imshow(panels{k}, []);
        colormap(ax, panelCmaps{k});
    end
    title(panelTitles{k}, 'Color','w', 'FontSize',8, 'Interpreter','none');
    set(ax, 'XColor','none', 'YColor','none');
end
sgtitle(figTitle, 'Color','w', 'FontSize',11, 'FontWeight','bold');
end
