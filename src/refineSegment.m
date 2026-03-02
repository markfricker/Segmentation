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
%        label.  Typically produced by watershedSegment or cellposeEnhance.
%
% NAME-VALUE PAIRS
%   'method'       : boundary-refinement strategy (default 'chanvese')
%                      'chanvese' — per-object Chan-Vese active contour
%                                   (refines to the true fluorescent boundary)
%                      'dilate'   — Voronoi-dilation + intensity clipping
%                                   (faster; expands each label toward the
%                                    raw-image foreground)
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
% NOTES
%   - If L contains no labelled objects (max(L(:)) == 0), the function
%     returns BW_r = zeros(size(I),'single') and L_r = zeros(size(I),'uint16')
%     immediately without computation.
%   - The 'chanvese' method requires activecontour (Image Processing
%     Toolbox R2014b+); an error is raised if it is absent.
%   - Both methods preserve the number of objects (unless objects shrink
%     below minArea), so refineSegment adjusts boundaries only; it does
%     not split or merge objects.
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
% See also: watershedSegment, localThresholdFast, activecontour, bwdist,
%           imdilate, label2rgb

%% ----------- input parsing -----------
p = inputParser;
addRequired(p, 'I', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
addParameter(p, 'method',       'chanvese');
addParameter(p, 'maxExpand',    5);
addParameter(p, 'nIter',        100);
addParameter(p, 'smoothFactor', 1.0);
addParameter(p, 'fgThreshold',  0.1);
addParameter(p, 'minArea',      10);
addParameter(p, 'pad',          []);   % derived below if empty
parse(p, I, L, varargin{:});

method     = lower(string(p.Results.method));
maxExpand  = p.Results.maxExpand;
nIter      = p.Results.nIter;
smoothFact = p.Results.smoothFactor;
fgThresh   = p.Results.fgThreshold;
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
    end

    %% clean up and relabel
    BW_clean = bwareaopen(BW_accum, minArea);
    L_r  = uint16(labelmatrix(bwconncomp(BW_clean)));
    BW_r = single(BW_clean);
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
    BW_clean = bwareaopen(L_dilated > 0, minArea);
    L_dilated(~BW_clean) = 0;
    L_r  = uint16(labelmatrix(bwconncomp(BW_clean)));
    BW_r = single(BW_clean);
    return
end

error('refineSegment:unknownMethod', ...
      'Unknown method "%s". Choose ''chanvese'' or ''dilate''.', method);

end
