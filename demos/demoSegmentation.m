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
% PIPELINE — Real mitochondrial image (Figs 4-6)
%   Load mitImage.mat from BlobFilters demos (set blobFiltersPath below).
%   Enhancement: rodGranulometryEnhance if BlobFilters is on path;
%                otherwise use the raw image directly.
%   localThresholdFast / watershedSegment applied to enhancement map.
%
% FIGURES PRODUCED
%   Fig 1 — Synthetic: thresholding  (Raw | Sauvola | Bernsen | Bradley | Niblack)
%   Fig 2 — Synthetic: watershed     (Raw | Sauvola mask | Dist | Grad | Marker)
%   Fig 3 — Synthetic: zoom on marker-WS
%   Fig 4 — Real:      thresholding  (enhancement map | Sauvola | Bernsen | Bradley | Niblack)
%   Fig 5 — Real:      watershed     (Raw | enhancement | Dist | Grad | Marker)
%   Fig 6 — Real:      zoom on marker-WS
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

% --- localThresholdFast (real mito image, objects ~8-12 px) ---------------
win_r = 31;

pSauv_r.windowsize = win_r; pSauv_r.k = 0.30; pSauv_r.r = 0.5;
pBern_r.windowsize = win_r; pBern_r.contrast = 0.05;
pBrad_r.windowsize = win_r; pBrad_r.k = 0.12;
pNibl_r.windowsize = win_r; pNibl_r.k = -0.15;

pWS_dist_r  = struct('method','distance', 'threshold',0.30, 'hMinima',2.0,  'minArea',50);
pWS_grad_r  = struct('method','gradient', 'threshold',0.30, 'smoothSigma',1.5, 'hMinima',0.04, 'minArea',50);
pWS_mark_r  = struct('method','marker',   'threshold',0.30, 'smoothSigma',1.5, 'hMinima',0.05, 'minArea',50);

% --- rodGranulometryEnhance (real image only) ----------------------------
pRod.lengths      = [8 12 16 20 28 36];
pRod.orientations = 8;
pRod.normalize    = true;

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
% Try mitImage.mat first (uint16 variable I), then Ireal.png.
mitoMatPath = fullfile(blobImagePath, 'mitImage.mat');
mitoImgPath = fullfile(blobImagePath, 'Ireal.png');

if ~isempty(blobFiltersPath) && exist(mitoMatPath, 'file')
    fprintf('Loading mitImage.mat...\n');
    tmp   = load(mitoMatPath, 'I');
    Ireal = im2single(tmp.I);
elseif ~isempty(blobFiltersPath) && exist(mitoImgPath, 'file')
    fprintf('Loading Ireal.png...\n');
    Ireal = im2single(imread(mitoImgPath));
    if size(Ireal,3) > 1
        Ireal = rgb2gray(Ireal);
    end
else
    fprintf('  Real image not found at:\n    %s\n    %s\n', mitoMatPath, mitoImgPath);
    fprintf('  Check blobImagePath at the top of this script.\n');
    fprintf('\nDone. 3 figures generated (synthetic only).\n');
    return
end

fprintf('  Size: %d x %d px\n', size(Ireal,2), size(Ireal,1));

% -------------------------------------------------------------------------
% B2.  Enhancement (rodGranulometryEnhance if available, else raw)
% -------------------------------------------------------------------------
if hasBlobFilters
    fprintf('  rodGranulometryEnhance... '); tic;
    Renh = rodGranulometryEnhance(Ireal, pRod);
    fprintf('%.2fs\n', toc);
    enhTitle = 'RodGran enhancement';
else
    fprintf('  BlobFilters not found — using raw image as enhancement map.\n');
    Renh     = Ireal;
    enhTitle = 'Raw (no enhancement)';
end

% -------------------------------------------------------------------------
% B3.  Thresholding on real enhancement map
% -------------------------------------------------------------------------
fprintf('\n--- Thresholding (real) ---\n');

fprintf('  sauvola...   '); tic;
BW_sauv_r = localThresholdFast(Renh, 'method','sauvola', ...
    'windowsize',pSauv_r.windowsize, 'k',pSauv_r.k, 'r',pSauv_r.r);
fprintf('%.2fs  (%d px)\n', toc, nnz(BW_sauv_r));

fprintf('  bernsen...   '); tic;
BW_bern_r = localThresholdFast(Renh, 'method','bernsen', ...
    'windowsize',pBern_r.windowsize, 'contrast',pBern_r.contrast);
fprintf('%.2fs  (%d px)\n', toc, nnz(BW_bern_r));

fprintf('  bradley...   '); tic;
BW_brad_r = localThresholdFast(Renh, 'method','bradley', ...
    'windowsize',pBrad_r.windowsize, 'k',pBrad_r.k);
fprintf('%.2fs  (%d px)\n', toc, nnz(BW_brad_r));

fprintf('  niblack...   '); tic;
BW_nibl_r = localThresholdFast(Renh, 'method','niblack', ...
    'windowsize',pNibl_r.windowsize, 'k',pNibl_r.k);
fprintf('%.2fs  (%d px)\n', toc, nnz(BW_nibl_r));

% -------------------------------------------------------------------------
% B4.  Watershed on real enhancement map
% -------------------------------------------------------------------------
fprintf('\n--- Watershed (real) ---\n');

fprintf('  distance...  '); tic;
[~, L_dist_r] = watershedSegment(Renh, 'method',pWS_dist_r.method, ...
    'threshold',pWS_dist_r.threshold, 'hMinima',pWS_dist_r.hMinima, ...
    'minArea',pWS_dist_r.minArea);
fprintf('%.2fs  (%d objects)\n', toc, max(L_dist_r(:)));

fprintf('  gradient...  '); tic;
[~, L_grad_r] = watershedSegment(Renh, 'method',pWS_grad_r.method, ...
    'threshold',pWS_grad_r.threshold, 'smoothSigma',pWS_grad_r.smoothSigma, ...
    'hMinima',pWS_grad_r.hMinima, 'minArea',pWS_grad_r.minArea);
fprintf('%.2fs  (%d objects)\n', toc, max(L_grad_r(:)));

fprintf('  marker...    '); tic;
[~, L_mark_r] = watershedSegment(Renh, 'method',pWS_mark_r.method, ...
    'threshold',pWS_mark_r.threshold, 'smoothSigma',pWS_mark_r.smoothSigma, ...
    'hMinima',pWS_mark_r.hMinima, 'minArea',pWS_mark_r.minArea);
fprintf('%.2fs  (%d objects)\n', toc, max(L_mark_r(:)));

% -------------------------------------------------------------------------
% B5.  Figures 4-6 — real image
% -------------------------------------------------------------------------
rgb_sauv_r = maskRGB(BW_sauv_r);
rgb_dist_r = labelRGB(L_dist_r);
rgb_grad_r = labelRGB(L_grad_r);
rgb_mark_r = labelRGB(L_mark_r);

figure(4);
set(gcf,'Name','Real — Thresholding','NumberTitle','off', ...
        'Color','k','Position',[30 680 1400 300]);
segFillFigure(4, {Renh, single(BW_sauv_r), single(BW_bern_r), ...
                  single(BW_brad_r), single(BW_nibl_r)}, ...
    sprintf('Thresholding methods — real mito (%s)', enhTitle), ...
    {enhTitle,'Sauvola','Bernsen','Bradley','Niblack'}, ...
    {'hot','gray','gray','gray','gray'});

figure(5);
set(gcf,'Name','Real — Watershed','NumberTitle','off', ...
        'Color','k','Position',[30 380 1400 300]);
segFillFigure(5, {Ireal, Renh, rgb_dist_r, rgb_grad_r, rgb_mark_r}, ...
    sprintf('Watershed methods — real mito (%s)', enhTitle), ...
    {'Raw', enhTitle, ...
     sprintf('Distance-WS (n=%d)', max(L_dist_r(:))), ...
     sprintf('Gradient-WS (n=%d)', max(L_grad_r(:))), ...
     sprintf('Marker-WS   (n=%d)', max(L_mark_r(:)))}, ...
    {'gray','hot',[],[],[]});

figure(6);
set(gcf,'Name','Real — Zoom marker-WS','NumberTitle','off', ...
        'Color','k','Position',[30 60 1400 320]);
segFillFigure(6, ...
    {imcrop(Ireal,roi_real1), imcrop(Renh,roi_real1), ...
     imcrop(rgb_sauv_r,roi_real1), imcrop(rgb_mark_r,roi_real1), ...
     imcrop(Ireal,roi_real2), imcrop(Renh,roi_real2), ...
     imcrop(rgb_sauv_r,roi_real2), imcrop(rgb_mark_r,roi_real2)}, ...
    sprintf('Zoom clusters — marker-WS  (n_total=%d)  |  roi1 (left) / roi2 (right)', ...
            max(L_mark_r(:))), ...
    {'Raw roi1',sprintf('%s roi1',enhTitle),'Sauvola roi1','Marker-WS roi1', ...
     'Raw roi2',sprintf('%s roi2',enhTitle),'Sauvola roi2','Marker-WS roi2'}, ...
    {'gray','hot',[],[],'gray','hot',[],[]});

fprintf('\nDone. 6 figures generated.\n');


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
