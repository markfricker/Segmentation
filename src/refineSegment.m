function [BW_r, L_r] = refineSegment(I, L, varargin)
% refineSegment  Boundary refinement of an instance segmentation
%
% USAGE
%   [BW_r, L_r] = refineSegment(I, L)
%   [BW_r, L_r] = refineSegment(I, L, Name, Value, ...)
%
% INPUTS
%   I  : 2-D raw fluorescence image (any numeric class; converted to
%        single [0,1] internally via im2single).  Must be the same size
%        as L.  Used as the intensity guide for boundary expansion.
%   L  : Instance label image, uint16, same size as I.
%        Background = 0; each object carries a unique positive integer
%        label.  Typically produced by watershedSegment or cellposeSegment.
%
% NAME-VALUE PAIRS
%   'method'       : boundary-refinement strategy (default 'chanvese')
%                      'chanvese' — per-object Chan-Vese active contour
%                                   (refines to the true fluorescent boundary)
%                      'dilate'   — Voronoi-dilation + intensity clipping
%                                   (faster; expands each label toward the
%                                    raw-image foreground)
%                      'halfmax'  — per-object adaptive threshold between
%                                   local background and peak
%                      'otsu'     — Voronoi-dilation + global Otsu threshold
%                      'grabcut'  — per-object GrabCut: the seed mask is the
%                                   foreground scribble, an outer band of the
%                                   dilated territory is the background
%                                   scribble, and the graph-cut solves the
%                                   true boundary in between
%   'maxExpand'    : maximum allowed expansion from the original object
%                    boundary, in pixels (default 5).
%                    Prevents adjacent objects from merging.
%   'nIter'        : number of active-contour iterations for the 'chanvese'
%                    method (default 100). Increase for larger objects or
%                    finer convergence; decrease for speed.
%   'smoothFactor' : Chan-Vese contour smoothness term weight (default 1.0).
%                    Higher values produce smoother contours; lower values
%                    allow finer boundary detail.
%   'fgThreshold'  : minimum raw intensity for pixel acceptance in the
%                    'dilate' method (default 0.1 on the [0,1] scale).
%                    Pixels below this level are treated as background even
%                    when they fall within the dilation radius.
%   'bgMargin'     : fraction of maxExpand (default 0.5) that stays an
%                    "undetermined" band next to the seed for the
%                    'grabcut' method.  Only the outer rim of the dilated
%                    territory beyond this margin is given to GrabCut as a
%                    hard background scribble.  Values that leave no
%                    usable gap (too close to 0, so the two scribbles
%                    would touch, or too close to 1, so the background
%                    scribble vanishes) fall back to the seed unchanged
%                    for that object -- GrabCut needs a real gap between
%                    its foreground and background evidence.
%   'minArea'      : minimum object area in pixels; objects smaller than
%                    this are discarded from the output (default 10).
%   'pad'          : bounding-box padding around each object patch used
%                    by the 'chanvese' method (default = maxExpand + 5).
%                    Must be >= maxExpand to avoid clipping the expansion
%                    zone.
%
% OUTPUTS
%   BW_r : refined binary mask, single [0,1], same size as I.
%          Foreground pixels = 1, background = 0.
%   L_r  : refined label image, uint16, same size as I.
%          Each object carries a unique positive integer; background = 0.
%          Object label IDs may differ from L after relabelling.
%          Use bwconncomp / regionprops on L_r directly.
%
% OVERVIEW
%   Watershed and classical-enhancer pipelines tend to produce tightly-
%   cropped binary masks that hug the peak of the enhancement response
%   rather than the true fluorescent boundary.  This leads to high
%   precision but low recall when benchmarked against Cellpose (which
%   labels the full organelle body).  refineSegment expands each object
%   mask outward from its watershed boundary to better capture the
%   complete organelle extent.
%
%   CHAN-VESE METHOD ('chanvese')  — per-object, intensity-guided
%     For each labelled instance k independently:
%       1. Extract a padded bounding-box patch from I.
%       2. Initialise the active contour as the input mask dilated by
%          maxExpand/2 pixels (giving the contour room to search outward
%          while starting close to the known boundary).
%       3. Run activecontour() with Chan-Vese energy (region statistics:
%          inside vs outside mean intensity) for nIter iterations.
%       4. Clip the refined contour to a maxExpand-radius dilation of the
%          original mask (hard spatial limit to prevent the contour from
%          jumping to distant bright regions or adjacent objects).
%       5. Accumulate refined masks; a pixel claimed by an earlier object
%          is not overwritten by a later one (prevents merging).
%     Requires: Image Processing Toolbox (activecontour, R2014b+).
%
%   DILATE METHOD ('dilate')  — global, Voronoi-partitioned
%     Single-pass algorithm; no per-object loop:
%       1. For every pixel, compute its distance to the nearest original
%          seed pixel and the label of that seed (bwdist Voronoi).
%       2. Accept an expanded pixel if BOTH:
%             distance to nearest seed <= maxExpand  (spatial limit)
%             I(pixel) >= fgThreshold               (intensity gate)
%          OR the pixel was already labelled in L    (originals preserved).
%       3. bwareaopen + relabel.
%     The Voronoi partition ensures adjacent objects never merge: each
%     candidate pixel is assigned to its nearest seed, so expansion stops
%     at the Voronoi boundary between two neighbouring objects.
%     Much faster than Chan-Vese (O(N) vs O(N*K) where K = object count).
%
%   GRABCUT METHOD ('grabcut')  — per-object, graph-cut with scribbles
%     For each labelled instance k independently:
%       1. Build the same exclusive Voronoi territory as 'halfmax' (capped
%          at maxExpand) — this is both the GrabCut region-of-interest and
%          the source of the background scribble.
%       2. Foreground scribble = the object's own seed mask (hard
%          constraint: these pixels are always foreground).
%       3. Background scribble = the part of the territory beyond a
%          maxExpand*bgMargin buffer around the seed (hard constraint:
%          always background). The band between the seed and this buffer
%          is left undetermined for the graph-cut to resolve using image
%          evidence, rather than being told outright.
%       4. superpixels() oversegments the padded patch; grabcut() solves
%          for the true object boundary within the territory, seeded by
%          the two scribbles.
%       5. Result is hard-clipped to the territory (prevents merging,
%          matching every other method here).
%     If the territory is too tight to leave any background band (i.e.
%     the buffer already covers the whole territory), the seed is
%     returned unchanged for that object — there is no background
%     evidence for GrabCut to learn from.
%     Requires: Image Processing Toolbox (grabcut, superpixels, R2018a+).
%     Slower than 'dilate'/'halfmax'/'otsu' (per-object superpixel +
%     graph-cut cost); similar cost profile to 'chanvese'.
%
% NOTES
%   - If L contains no labelled objects (max(L(:)) == 0), the function
%     returns BW_r = zeros(size(I),'single') and L_r = zeros(size(I),'uint16')
%     immediately without computation.
%   - The 'chanvese' method requires activecontour (Image Processing
%     Toolbox R2014b+); an error is raised if it is absent.
%   - The 'grabcut' method requires grabcut/superpixels (Image Processing
%     Toolbox R2018a+); an error is raised if grabcut is absent.
%   - All methods preserve object identity: adjacent objects are never
%     merged.  Chan-Vese and GrabCut use a no-overwrite accumulator and
%     per-label area filtering; the Voronoi-based methods (dilate, otsu,
%     halfmax, and GrabCut's territory) use the disjoint partition
%     directly with per-label area filtering.  No method routes the final
%     label image through bwconncomp (which would merge touching objects).
%   - Setting maxExpand = 0 with 'dilate' returns exactly the original
%     foreground (no expansion).  With 'chanvese' it still runs the AC
%     (which may contract slightly) but clips the result to the original
%     mask.
%   - For very elongated mitochondria, 'dilate' with maxExpand = 5-8 and
%     fgThreshold = 0.08-0.12 typically gives the best speed/quality
%     trade-off.  Use 'chanvese' when precise boundary localisation
%     matters more than runtime.
%
% REFERENCES
%   Chan, T.F. and Vese, L.A. (2001). Active contours without edges.
%   IEEE Transactions on Image Processing, 10(2), 266-277.
%
%   Soille, P. (2003). Morphological Image Analysis: Principles and
%   Applications. 2nd ed. Springer.  (Voronoi dilation, pp. 81-83.)
%
% EXAMPLE
%   % Two-stage pipeline: marker watershed → Chan-Vese refinement
%   R = capsuleEnhance(Ism, pCap);
%   [~, L_ws] = watershedSegment(R, 'threshold', 0.3, 'hMinima', 0.05);
%   [BW_r, L_r] = refineSegment(Iraw, L_ws, ...
%                                'method',    'chanvese', ...
%                                'maxExpand', 6, ...
%                                'nIter',     150);
%   figure;
%   subplot(1,3,1); imshow(Iraw, []);                  title('Raw');
%   subplot(1,3,2); imshow(label2rgb(L_ws,'jet','k')); title('Watershed');
%   subplot(1,3,3); imshow(label2rgb(L_r, 'jet','k')); title('Refined (CV)');
%
%   % Fast Voronoi-dilation refinement
%   [BW_r, L_r] = refineSegment(Iraw, L_ws, ...
%                                'method',       'dilate', ...
%                                'maxExpand',    8, ...
%                                'fgThreshold',  0.10);
%
%   % GrabCut refinement (seed = foreground scribble, outer dilation
%   % ring = background scribble)
%   [BW_r, L_r] = refineSegment(Iraw, L_ws, ...
%                                'method',    'grabcut', ...
%                                'maxExpand', 8, ...
%                                'bgMargin',  0.5);
%
% See also: watershedSegment, cellposeSegment, localThresholdFast,
%           activecontour, grabcut, superpixels, bwdist, imdilate, label2rgb

%% ----------- input parsing -----------
p = inputParser;
addRequired(p, 'I', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
addParameter(p, 'method',       'chanvese');
addParameter(p, 'maxExpand',    5);
addParameter(p, 'nIter',        100);
addParameter(p, 'smoothFactor', 1.0);
addParameter(p, 'fgThreshold',  0.1);
addParameter(p, 'level',        0.5);  % 'halfmax': fraction between local bg and peak
addParameter(p, 'bgMargin',     0.5);  % 'grabcut': fraction of maxExpand kept undetermined next to the seed
addParameter(p, 'minArea',      10);
addParameter(p, 'pad',          []);   % derived below if empty
parse(p, I, L, varargin{:});

method     = lower(string(p.Results.method));
maxExpand  = p.Results.maxExpand;
nIter      = p.Results.nIter;
smoothFact = p.Results.smoothFactor;
fgThresh   = p.Results.fgThreshold;
level      = p.Results.level;
bgMargin   = p.Results.bgMargin;
minArea    = p.Results.minArea;
padSz      = p.Results.pad;

if isempty(padSz)
    padSz = maxExpand + 5;
end

I = im2single(I);
L = uint16(L);

[imgH, imgW] = size(I);

%% fast exit for empty label image
if max(L(:)) == 0
    BW_r = zeros(imgH, imgW, 'single');
    L_r  = zeros(imgH, imgW, 'uint16');
    return
end

%% ============================================================
%% Chan-Vese per-object active contour
%% ============================================================
if method == "chanvese"

    if ~exist('activecontour', 'file')
        error('refineSegment:noActiveContour', ...
              ['activecontour() not found.  The ''chanvese'' method ' ...
               'requires the Image Processing Toolbox (R2014b+).']);
    end

    nLabels  = double(max(L(:)));
    BW_accum = false(imgH, imgW);   % accumulator — earlier objects take priority
    L_accum  = zeros(imgH, imgW, 'uint16');  % per-pixel label accumulator

    % structuring elements reused across objects
    seInit = strel('disk', max(1, round(maxExpand / 2)));
    seClip = strel('disk', max(0, maxExpand));

    for k = 1:nLabels
        mask_k = (L == k);
        if ~any(mask_k(:))
            continue   % label may be absent (gap after upstream filtering)
        end

        % --- padded bounding box (clamped to image) ---
        [rows, cols] = find(mask_k);
        r1 = max(1,    min(rows) - padSz);
        r2 = min(imgH, max(rows) + padSz);
        c1 = max(1,    min(cols) - padSz);
        c2 = min(imgW, max(cols) + padSz);

        I_patch    = I(r1:r2, c1:c2);
        mask_patch = mask_k(r1:r2, c1:c2);

        % --- initial contour: dilated mask gives AC room to search outward ---
        init_patch = imdilate(mask_patch, seInit);

        % --- hard clipping envelope: original mask dilated by maxExpand ---
        clip_patch = imdilate(mask_patch, seClip);

        % --- Chan-Vese active contour (region-based, inside vs outside mean) ---
        try
            refined_patch = activecontour(I_patch, init_patch, nIter, ...
                                          'Chan-Vese', ...
                                          'SmoothFactor', smoothFact);
        catch
            % Fallback: use dilated initialisation as-is on AC failure
            refined_patch = init_patch;
        end

        % --- apply spatial clip ---
        refined_patch = refined_patch & clip_patch;

        % --- write back: earlier objects take priority, no overwrites ---
        row_range = r1:r2;
        col_range = c1:c2;
        new_pixels = refined_patch & ~BW_accum(row_range, col_range);
        BW_accum(row_range, col_range) = ...
            BW_accum(row_range, col_range) | new_pixels;
        % Record which label owns each new pixel (preserves boundaries)
        L_patch = L_accum(row_range, col_range);
        L_patch(new_pixels) = uint16(k);
        L_accum(row_range, col_range) = L_patch;
    end

    %% clean up and relabel
    % Remove per-object labels below minArea (per-label check avoids
    % bwconncomp, which would merge touching objects into one component).
    for k = 1:nLabels
        if nnz(L_accum == k) < minArea
            L_accum(L_accum == k) = 0;
        end
    end
    uL = setdiff(unique(L_accum(:)), uint16(0));
    L_r = zeros(imgH, imgW, 'uint16');
    for ki = 1:numel(uL)
        L_r(L_accum == uL(ki)) = ki;
    end
    BW_r = single(L_r > 0);
    return
end

%% ============================================================
%% Voronoi-dilation + intensity gating  (fast global method)
%% ============================================================
if method == "dilate"

    seedMask = L > 0;

    % For each pixel: distance to nearest seed and its linear index
    [D, nearest_idx] = bwdist(seedMask);

    % Propagate labels via nearest-seed assignment (Voronoi partition)
    L_voronoi = uint16(L(nearest_idx));

    % Accept expanded pixels that are:
    %   (a) within the expansion radius AND above the intensity gate, OR
    %   (b) already labelled in the original mask (never shrink existing labels)
    valid    = (D <= maxExpand & I >= fgThresh) | seedMask;
    L_dilated = L_voronoi;
    L_dilated(~valid) = 0;

    %% clean up and relabel
    % Use per-label area filtering rather than bwconncomp, which would
    % merge adjacent Voronoi cells that share a boundary edge.
    uL = setdiff(unique(L_dilated(:)), uint16(0));
    L_r    = zeros(imgH, imgW, 'uint16');
    ki_out = uint16(0);
    for ki = 1:numel(uL)
        mask_k = (L_dilated == uL(ki));
        if nnz(mask_k) >= minArea
            ki_out = ki_out + 1;
            L_r(mask_k) = ki_out;
        end
    end
    BW_r = single(L_r > 0);
    return
end

%% ============================================================
%% Local half-max threshold within a non-merging Voronoi territory
%% ============================================================
if method == "halfmax"
    % Each object gets an exclusive expanded territory (Voronoi, capped at
    % maxExpand) and is then thresholded at LEVEL of the way from its local
    % background to its peak — an intensity-adaptive boundary per object,
    % rather than one fixed threshold for all.  Adjacent objects cannot merge
    % (territories are disjoint).
    seedMask = L > 0;
    [D, nearest_idx] = bwdist(seedMask);
    L_terr = uint16(L(nearest_idx));
    L_terr(D > maxExpand) = 0;                  % exclusive territory per label

    nLabels = double(max(L(:)));
    L_accum = zeros(imgH, imgW, 'uint16');

    for k = 1:nLabels
        objMask  = (L == k);
        terrMask = (L_terr == k);
        if ~any(objMask(:)) || ~any(terrMask(:)), continue; end

        peak = pctl(I(objMask), 95);
        ring = terrMask & ~objMask;             % expansion zone ~ local background
        if any(ring(:))
            bg = pctl(I(ring), 50);
        else
            bg = pctl(I(terrMask), 10);
        end
        thr = bg + level * (peak - bg);

        accept = (terrMask & I >= thr);
        seed   = objMask & accept;
        if ~any(seed(:)), seed = objMask; end
        refined = imreconstruct(seed, accept | objMask);   % keep core, grow to half-max

        new = refined & (L_accum == 0);         % no-overwrite (territories disjoint)
        L_accum(new) = uint16(k);
    end

    for k = 1:nLabels
        if nnz(L_accum == k) < minArea
            L_accum(L_accum == k) = 0;
        end
    end
    uL  = setdiff(unique(L_accum(:)), uint16(0));
    L_r = zeros(imgH, imgW, 'uint16');
    for ki = 1:numel(uL)
        L_r(L_accum == uL(ki)) = ki;
    end
    BW_r = single(L_r > 0);
    return
end

%% ============================================================
%% Voronoi-dilation + global Otsu threshold  (automatic intensity gate)
%% ============================================================
if method == "otsu"
    % Compute global Otsu threshold from the intensity image.
    % graythresh expects a normalised [0,1] image (which im2single gives us).
    thr = double(graythresh(I));

    seedMask = L > 0;
    [D, nearest_idx] = bwdist(seedMask);
    L_voronoi = uint16(L(nearest_idx));

    valid     = (D <= maxExpand & I >= thr) | seedMask;
    L_dilated = L_voronoi;
    L_dilated(~valid) = 0;

    uL = setdiff(unique(L_dilated(:)), uint16(0));
    L_r    = zeros(imgH, imgW, 'uint16');
    ki_out = uint16(0);
    for ki = 1:numel(uL)
        mask_k = (L_dilated == uL(ki));
        if nnz(mask_k) >= minArea
            ki_out = ki_out + 1;
            L_r(mask_k) = ki_out;
        end
    end
    BW_r = single(L_r > 0);
    return
end

%% ============================================================
%% GrabCut with seed=foreground / dilation-ring=background scribbles
%% ============================================================
if method == "grabcut"

    if ~exist('grabcut', 'file')
        error('refineSegment:noGrabcut', ...
              ['grabcut() not found.  The ''grabcut'' method requires ' ...
               'the Image Processing Toolbox (R2018a+).']);
    end

    % Exclusive Voronoi territory per label, capped at maxExpand -- same
    % construction as 'halfmax'.  This is both the GrabCut region of
    % interest (the zone GrabCut is allowed to grow into) and the source
    % of the background scribble.
    seedMask = L > 0;
    [D, nearest_idx] = bwdist(seedMask);
    L_terr = uint16(L(nearest_idx));
    L_terr(D > maxExpand) = 0;

    nLabels = double(max(L(:)));
    L_accum = zeros(imgH, imgW, 'uint16');

    % bufRadius < 1 leaves no gap between the foreground and background
    % scribbles; since GrabCut's graph-cut operates on whole superpixels
    % (not individual pixels), a superpixel straddling that zero-width
    % boundary would receive contradictory hard constraints and can drop
    % genuine seed pixels rather than simply refuse to expand. Treat it
    % the same as "no background evidence": skip GrabCut, keep the seed.
    bufRadius = round(maxExpand * bgMargin);
    hasBuffer = bufRadius >= 1;
    if hasBuffer
        seBuffer = strel('disk', bufRadius);
    end

    for k = 1:nLabels
        objMask  = (L == k);
        terrMask = (L_terr == k);
        if ~any(objMask(:)) || ~any(terrMask(:)), continue; end

        % --- padded bounding box around the territory (clamped) ---
        [rows, cols] = find(terrMask);
        r1 = max(1,    min(rows) - 2);
        r2 = min(imgH, max(rows) + 2);
        c1 = max(1,    min(cols) - 2);
        c2 = min(imgW, max(cols) + 2);

        I_patch    = I(r1:r2, c1:c2);
        obj_patch  = objMask(r1:r2, c1:c2);
        terr_patch = terrMask(r1:r2, c1:c2);

        bg_patch = false(size(obj_patch));
        if hasBuffer
            % background scribble: the outer rim of the territory beyond
            % a bufRadius buffer around the seed.  The band between the
            % seed and this buffer is left undetermined for GrabCut to
            % resolve from image evidence rather than being told outright.
            buf_patch = imdilate(obj_patch, seBuffer);
            bg_patch  = terr_patch & ~buf_patch;
        end

        if ~any(bg_patch(:))
            % No room for a background scribble (bgMargin too high, or
            % too low to leave a safe gap, for this territory) -- nothing
            % for GrabCut to learn from.
            refined_patch = obj_patch;
        else
            nSP = max(20, round(numel(I_patch) / 9));
            Lsp = superpixels(I_patch, nSP);
            refined_patch = grabcut(I_patch, Lsp, terr_patch, ...
                find(obj_patch), find(bg_patch));
            refined_patch = refined_patch & terr_patch;   % hard clip (belt-and-braces)
        end

        % --- write back: earlier objects take priority, no overwrites ---
        row_range = r1:r2;
        col_range = c1:c2;
        new = refined_patch & (L_accum(row_range, col_range) == 0);
        L_patch = L_accum(row_range, col_range);
        L_patch(new) = uint16(k);
        L_accum(row_range, col_range) = L_patch;
    end

    for k = 1:nLabels
        if nnz(L_accum == k) < minArea
            L_accum(L_accum == k) = 0;
        end
    end
    uL  = setdiff(unique(L_accum(:)), uint16(0));
    L_r = zeros(imgH, imgW, 'uint16');
    for ki = 1:numel(uL)
        L_r(L_accum == uL(ki)) = ki;
    end
    BW_r = single(L_r > 0);
    return
end

error('refineSegment:unknownMethod', ...
      'Unknown method "%s". Choose ''chanvese'', ''dilate'', ''halfmax'', ''otsu'' or ''grabcut''.', method);

end

% -------------------------------------------------------------------------
function v = pctl(x, q)
%PCTL  Percentile without the Statistics Toolbox (q in 0..100).
    if isempty(x), v = 0; return; end
    xs = sort(x(:));
    v  = xs(min(numel(xs), max(1, round(q/100 * numel(xs)))));
end
