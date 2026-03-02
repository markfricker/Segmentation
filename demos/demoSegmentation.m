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
% PIPELINE — Real mitochondrial image (Figs 4-7)
%   Load mitImage.mat from BlobFilters demos (set blobImagePath below).
%   Preprocessing : orientedGaussSmooth  → Ism  (if BlobFilters available)
%   Enhancers (on Ism): logEnhance | fiberEnhance | capsuleEnhance(DoC)
%                       | rodGranulometryEnhance
%   Cellpose (on raw):  cellposeEnhance  → binary mask + label image
%   Segmentation: marker-controlled watershedSegment on each enhancement map.
%   Comparison: watershed labels from each enhancer vs Cellpose labels.
%
% FIGURES PRODUCED
%   Fig 1 — Synthetic: thresholding  (Raw | Sauvola | Bernsen | Bradley | Niblack)
%   Fig 2 — Synthetic: watershed     (Raw | Sauvola mask | Dist | Grad | Marker)
%   Fig 3 — Synthetic: zoom on marker-WS
%   Fig 4 — Real:      enhancement maps  (Raw | log | fib | cap | rod | CP mask | CP labels)
%   Fig 5 — Real:      watershed vs CP   (Raw | log-WS | fib-WS | cap-WS | rod-WS | CP labels)
%   Fig 6 — Real:      zoom cluster 1    (same 6 columns as Fig 5)
%   Fig 7 — Real:      zoom cluster 2    (same 6 columns as Fig 5)
%
% REQUIREMENTS
%   localThresholdFast.m and watershedSegment.m (added automatically).
%   Image Processing Toolbox (bwdist, watershed, imextendedmax, etc.)
%   BlobFilters path for enhancement (optional; raw image used as fallback).

clear; clc; close all;

% =========================================================================
% CONFIGURATION  ← adjust these paths for your machine
% =========================================================================
% Directory containing BlobFilters enhancer .m files (rodGranulometryEnhance etc.)
% The main repo has these at the repo root (not in a src/ subdirectory).
% Set to '' to skip enhancement and use the raw image directly.
blobFunctionPath = 'C:\Users\dops0035\Documents\Research\Matlab Projects\BlobFilters_sandbox\src';

% Directory containing mitImage.mat and Ireal.png.
% These live in the demos/ folder of the active worktree.
blobImagePath = 'C:\Users\dops0035\Documents\Research\Matlab Projects\BlobFilters_sandbox\src\.claude\worktrees\thirsty-wescoff\demos';

% Zoom ROIs for real image:  [x_left, y_top, width, height]  (imcrop convention)
roi_real1 = [100  80  220 220];   % cluster 1  (adjust after Fig 5)
roi_real2 = [330 220  220 220];   % cluster 2  (adjust after Fig 5)

% Synthetic image zoom ROI
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

% --- watershedSegment (real mito, marker method — applied to all enhancers)
% A single consistent parameter set keeps the comparison fair.
pWS_r = struct('method','marker', 'threshold',0.30, ...
               'smoothSigma',1.5, 'hMinima',0.05, 'minArea',50);

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

% =========================================================================
% =========================================================================
%  PART B — REAL MITOCHONDRIAL IMAGE
% =========================================================================
% =========================================================================

fprintf('\n=== PART B: Real mitochondrial image ===\n');

% -------------------------------------------------------------------------
% B1.  Load real image
% -------------------------------------------------------------------------
mitoMatPath = fullfile(blobImagePath, 'mitImage.mat');
mitoImgPath = fullfile(blobImagePath, 'Ireal.png');

if ~isempty(blobImagePath) && exist(mitoMatPath, 'file')
    fprintf('Loading mitImage.mat...\n');
    tmp   = load(mitoMatPath, 'I');
    Ireal = im2single(tmp.I);
elseif ~isempty(blobImagePath) && exist(mitoImgPath, 'file')
    fprintf('Loading Ireal.png...\n');
    Ireal = im2single(imread(mitoImgPath));
    if size(Ireal,3) > 1, Ireal = rgb2gray(Ireal); end
else
    fprintf('  Real image not found at:\n    %s\n    %s\n', mitoMatPath, mitoImgPath);
    fprintf('  Check blobImagePath at the top of this script.\n');
    fprintf('\nDone. 3 figures generated (synthetic only).\n');
    return
end
fprintf('  Size: %d x %d px\n', size(Ireal,2), size(Ireal,1));

% -------------------------------------------------------------------------
% B2.  BlobFilters enhancer parameters  (same as demoMitoEnhance)
% -------------------------------------------------------------------------
pOGS.sigmaAlong   = 4;    pOGS.sigmaAcross  = 1.5;
pOGS.orientations = 8;    pOGS.sigmaGrad    = 1.5;  pOGS.sigmaInt = 5;

pLog.sigmas    = [2 3 4 5 6];   pLog.normalize = true;

pFib.widths    = [6 7 8 9 10];  pFib.multimode = 'stack';  pFib.normalize = true;

pCap.lengths      = [12 16 20 28 36 40];  pCap.width      = 8;
pCap.wideWidth    = 18;                    pCap.alpha      = 0.55;
pCap.orientations = 12;                    pCap.mode       = 'doc';
pCap.normalize    = true;

pRod.lengths      = [8 12 16 20 28 36];   pRod.orientations = 8;
pRod.normalize    = true;

pCP.model         = 'cyto3';   pCP.diameter      = 10;
pCP.cellProb      = -2;        pCP.flowThreshold = 0.8;

% -------------------------------------------------------------------------
% B3.  Preprocessing: OGS smooth
% -------------------------------------------------------------------------
if hasBlobFilters && exist('orientedGaussSmooth','file')
    fprintf('  orientedGaussSmooth... '); tic;
    Ism = orientedGaussSmooth(Ireal, pOGS);
    fprintf('%.2fs\n', toc);
else
    Ism = Ireal;   % fallback: no OGS
end

% -------------------------------------------------------------------------
% B4.  Classical enhancers (on Ism)
% -------------------------------------------------------------------------
fprintf('\n--- Classical enhancers (on OGS-smoothed image) ---\n');

if hasBlobFilters
    fprintf('  logEnhance...             '); tic;
    logR = logEnhance(Ism, pLog);
    fprintf('%.2fs\n', toc);

    fprintf('  fiberEnhance...           ');
    try
        tic; fibR = fiberEnhance(Ism, pFib); fprintf('%.2fs\n', toc);
    catch ME
        fprintf('SKIPPED (%s)\n', ME.message);
        fibR = zeros(size(Ism), 'single');
    end

    fprintf('  capsuleEnhance (DoC)...   '); tic;
    capR = capsuleEnhance(Ism, pCap);
    fprintf('%.2fs\n', toc);

    fprintf('  rodGranulometryEnhance... '); tic;
    rodR = rodGranulometryEnhance(Ism, pRod);
    fprintf('%.2fs\n', toc);
else
    fprintf('  BlobFilters not found — using raw image for all enhancer slots.\n');
    logR = Ireal;  fibR = Ireal;  capR = Ireal;  rodR = Ireal;
end

% -------------------------------------------------------------------------
% B5.  Cellpose (on raw image)
% -------------------------------------------------------------------------
fprintf('\n--- Cellpose ---\n');
hasCellpose = exist('cellpose','file') ~= 0;
if hasCellpose
    fprintf('  cellposeEnhance... ');
    try
        tic;
        [cpMask, cpL] = cellposeEnhance(Ireal, pCP);
        fprintf('%.2fs  (%d objects)\n', toc, max(cpL(:)));
    catch ME
        fprintf('SKIPPED (%s)\n', ME.message);
        cpMask = zeros(size(Ireal),'single');
        cpL    = zeros(size(Ireal),'uint16');
        hasCellpose = false;
    end
else
    fprintf('  SKIPPED (Cellpose add-on not installed)\n');
    cpMask = zeros(size(Ireal),'single');
    cpL    = zeros(size(Ireal),'uint16');
end

% -------------------------------------------------------------------------
% B6.  Marker watershed on each enhancement map
%      Single consistent parameter set — keeps the comparison fair.
% -------------------------------------------------------------------------
fprintf('\n--- Marker watershed on each enhancement map ---\n');

runWS = @(R) watershedSegment(R, 'method',pWS_r.method, ...
    'threshold',pWS_r.threshold, 'smoothSigma',pWS_r.smoothSigma, ...
    'hMinima',pWS_r.hMinima, 'minArea',pWS_r.minArea);

fprintf('  log-WS...   '); tic; [~, L_log] = runWS(logR); fprintf('%.2fs  (%d objects)\n', toc, max(L_log(:)));
fprintf('  fib-WS...   '); tic; [~, L_fib] = runWS(fibR); fprintf('%.2fs  (%d objects)\n', toc, max(L_fib(:)));
fprintf('  cap-WS...   '); tic; [~, L_cap] = runWS(capR); fprintf('%.2fs  (%d objects)\n', toc, max(L_cap(:)));
fprintf('  rod-WS...   '); tic; [~, L_rod] = runWS(rodR); fprintf('%.2fs  (%d objects)\n', toc, max(L_rod(:)));

% -------------------------------------------------------------------------
% B7.  Colour overlays
% -------------------------------------------------------------------------
rgb_log = labelRGB(L_log);
rgb_fib = labelRGB(L_fib);
rgb_cap = labelRGB(L_cap);
rgb_rod = labelRGB(L_rod);
rgb_cpL = labelRGB(cpL);

n_cp  = max(cpL(:));
n_log = max(L_log(:));  n_fib = max(L_fib(:));
n_cap = max(L_cap(:));  n_rod = max(L_rod(:));

% Shared panel titles / colormaps for Figs 5-7
wsTitles = {'Raw', ...
    sprintf('log-WS  (n=%d)',  n_log), ...
    sprintf('fib-WS  (n=%d)',  n_fib), ...
    sprintf('cap-WS  (n=%d)',  n_cap), ...
    sprintf('rod-WS  (n=%d)',  n_rod), ...
    sprintf('Cellpose (n=%d)', n_cp)};
wsCmaps  = {'gray',[],[],[],[],[]};

% -------------------------------------------------------------------------
% B8.  Figures 4-7
% -------------------------------------------------------------------------

% Fig 4 — enhancement maps
figure(4);
set(gcf,'Name','Real — Enhancement maps','NumberTitle','off', ...
        'Color','k','Position',[30 680 1820 300]);
segFillFigure(4, {Ireal, logR, fibR, capR, rodR, cpMask, rgb_cpL}, ...
    'Enhancement maps — real mitochondria', ...
    {'Raw','log','fib','cap (DoC)','rodGran', ...
     sprintf('CP mask (n=%d)',n_cp), sprintf('CP labels (n=%d)',n_cp)}, ...
    {'gray','hot','hot','hot','hot','hot',[]});

% Fig 5 — watershed vs Cellpose, full FOV
figure(5);
set(gcf,'Name','Real — Watershed vs Cellpose (full FOV)','NumberTitle','off', ...
        'Color','k','Position',[30 380 1820 300]);
segFillFigure(5, {Ireal, rgb_log, rgb_fib, rgb_cap, rgb_rod, rgb_cpL}, ...
    'Marker watershed on each enhancer vs Cellpose — full FOV', ...
    wsTitles, wsCmaps);

% Fig 6 — zoom cluster 1
figure(6);
set(gcf,'Name','Real — Zoom cluster 1','NumberTitle','off', ...
        'Color','k','Position',[30 60 1820 300]);
segFillFigure(6, ...
    {imcrop(Ireal,roi_real1), imcrop(rgb_log,roi_real1), imcrop(rgb_fib,roi_real1), ...
     imcrop(rgb_cap,roi_real1), imcrop(rgb_rod,roi_real1), imcrop(rgb_cpL,roi_real1)}, ...
    sprintf('Zoom cluster 1  [x=%d y=%d %dx%d px]', roi_real1(1),roi_real1(2),roi_real1(3),roi_real1(4)), ...
    wsTitles, wsCmaps);

% Fig 7 — zoom cluster 2
figure(7);
set(gcf,'Name','Real — Zoom cluster 2','NumberTitle','off', ...
        'Color','k','Position',[30 0 1820 300]);
segFillFigure(7, ...
    {imcrop(Ireal,roi_real2), imcrop(rgb_log,roi_real2), imcrop(rgb_fib,roi_real2), ...
     imcrop(rgb_cap,roi_real2), imcrop(rgb_rod,roi_real2), imcrop(rgb_cpL,roi_real2)}, ...
    sprintf('Zoom cluster 2  [x=%d y=%d %dx%d px]', roi_real2(1),roi_real2(2),roi_real2(3),roi_real2(4)), ...
    wsTitles, wsCmaps);

fprintf('\nDone. 7 figures generated.\n');


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
    title(panelTitles{k}, 'Color','w', 'FontSize',9);
    set(ax, 'XColor','none', 'YColor','none');
end
sgtitle(figTitle, 'Color','w', 'FontSize',12, 'FontWeight','bold');
end
