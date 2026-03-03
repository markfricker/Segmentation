function [BW, L] = watershedSegment(I, varargin)
% watershedSegment  Watershed-based segmentation with three flooding strategies
%
% USAGE
%   [BW, L] = watershedSegment(I)
%   [BW, L] = watershedSegment(I, Name, Value, ...)
%
% INPUTS
%   I  : 2-D grayscale or enhancement image (any numeric class; converted
%        to single [0,1] internally via im2single).
%        For the 'distance' method, I is binarised at 'threshold' first;
%        for 'gradient' and 'marker' methods, I is used as a continuous map.
%
% NAME-VALUE PAIRS
%   'method'      : watershed flooding strategy (default 'marker')
%                     'distance' — distance-transform watershed
%                     'gradient' — gradient-magnitude watershed
%                     'marker'   — marker-controlled watershed (default)
%   'threshold'   : foreground binarisation level in [0,1] (default 0.5)
%                   Acts as the HIGH threshold when 'threshLow' is also set.
%   'threshLow'   : lower threshold for two-pass hysteresis binarisation
%                   (default [], disabled).  When set, only the initial
%                   foreground mask is affected: pixels above 'threshold'
%                   unconditionally seed the foreground; pixels between
%                   'threshLow' and 'threshold' are included only if they
%                   belong to the same connected component as at least one
%                   high-threshold pixel.  Useful on enhancement maps with
%                   gradual boundaries (fiberEnhance, rodGranulometryEnhance)
%                   where a single low threshold over-segments and a single
%                   high threshold clips object edges.
%   'smoothSigma' : sigma of Gaussian pre-smoothing applied before
%                   gradient or distance computation (default 1.0)
%   'hMinima'     : h-minima suppression depth; controls over-segmentation
%                   by merging basins shallower than this value
%                   (default 0.02 for gradient/marker, 1.0 for distance)
%   'minArea'     : minimum object area in pixels; objects smaller than
%                   this are removed from the output (default 10)
%   'seedErosion' : radius of disk used to erode the foreground mask to
%                   generate foreground seed markers in the 'marker'
%                   method (default 2). Increase for larger objects.
%   'bgDilation'  : radius of disk used to dilate the foreground mask
%                   before inverting to obtain the background marker in
%                   the 'marker' method (default 3).
%
% OUTPUTS
%   BW : binary segmentation mask, single precision [0,1], same size as I.
%        Foreground pixels = 1, background = 0.  Consistent with the
%        output convention of cellposeEnhance and other BlobFilters tools.
%   L  : label image, uint16, same size as I.
%        Each detected object carries a unique positive integer label.
%        Background = 0.  Use label2rgb(L) for display.
%
% OVERVIEW
%   The watershed transform treats the image as a topographic relief map
%   and simulates flooding from local minima upward.  When two flood
%   basins meet, a watershed ridge (object boundary) is formed.  Applied
%   to a suitable relief image, each basin corresponds to one object.
%
%   The key challenge is over-segmentation: noise and texture create many
%   spurious minima, each spawning its own basin.  The three methods
%   address this differently:
%
%   DISTANCE TRANSFORM ('distance')
%     Steps:
%       1. Binarise I at 'threshold' to get foreground mask BW0.
%       2. Compute the Euclidean distance transform D = bwdist(~BW0).
%       3. Suppress shallow minima in -D with imhmin (depth 'hMinima').
%       4. Watershed on the suppressed -D: each local maximum of D
%          becomes a basin, corresponding to one object.
%       5. Constrain labels to foreground.
%     Best for: compact, roughly circular objects that are well separated
%     after thresholding.  Poor for elongated tubular structures where
%     the distance ridgeline is thin and fragmented.
%
%   GRADIENT MAGNITUDE ('gradient')
%     Steps:
%       1. Smooth I with a Gaussian (sigma = 'smoothSigma').
%       2. Compute gradient magnitude G = |grad(I_smooth)|.
%       3. Suppress shallow minima in G with imhmin.
%       4. Watershed on G: gradient valleys (object interiors) become
%          basins; gradient ridges (boundaries) become watershed lines.
%       5. Optionally mask result to foreground defined by 'threshold'.
%     Best for: images where object edges are sharper than their interiors,
%     such as strong membrane staining.  Tends to over-segment if the
%     gradient is not clean; increase 'hMinima' to compensate.
%
%   MARKER-CONTROLLED ('marker') — DEFAULT
%     Steps:
%       1. Binarise I at 'threshold' to get foreground mask BW0.
%       2. Foreground markers: imextendedmax on smoothed I, restricted to
%          BW0.  Each extended maximum becomes one seed (one future object).
%          'hMinima' controls how prominent a maximum must be to qualify.
%       3. Background marker: ~imdilate(BW0, disk('bgDilation')).
%       4. Build relief: negate smoothed I so object centres are deep.
%       5. Impose markers as forced minima on the relief (imimposemin).
%       6. Watershed: exactly one basin per foreground seed + background.
%       7. Constrain labels to foreground mask.
%     Best for: general purpose use on enhancement maps.  The marker step
%     decouples seed detection from boundary detection, giving explicit
%     control over the number of objects.  Recommended starting method.
%
%   POST-PROCESSING
%     All methods remove objects smaller than 'minArea' pixels using
%     bwareaopen before constructing the output label image.
%
% NOTES
%   - The watershed ridge pixels are assigned label 0 by MATLAB's
%     watershed() function.  These are absorbed into the nearest labelled
%     region during the final masking step so that BW has no gaps.
%   - For very thin tubular mitochondria, the 'marker' method with a small
%     'hMinima' (e.g. 0.01) and 'seedErosion' = 1 is recommended.
%   - For post-Cellpose use, prefer cellposeEnhance directly; this function
%     is intended for classical (non-deep-learning) pipelines.
%   - All three methods accept any [0,1]-normalised enhancement map as
%     input, including the outputs of fiberEnhance, rodGranulometryEnhance,
%     structureTensorEnhance, etc.
%
% REFERENCES
%   Beucher, S. and Lantuejoul, C. (1979). Use of watersheds in contour
%   detection. Proceedings of the International Workshop on Image
%   Processing, Real-Time Edge and Motion Detection.
%
%   Vincent, L. and Soille, P. (1991). Watersheds in digital spaces: an
%   efficient algorithm based on immersion simulations. IEEE Transactions
%   on Pattern Analysis and Machine Intelligence, 13(6), 583-598.
%
%   Beucher, S. and Meyer, F. (1993). The morphological approach to
%   segmentation: the watershed transformation. In: Mathematical
%   Morphology in Image Processing (ed. Dougherty, E.R.), pp. 433-481.
%   Marcel Dekker.
%
%   Soille, P. (2003). Morphological Image Analysis: Principles and
%   Applications. 2nd ed. Springer.
%
% EXAMPLE
%   % Marker-controlled watershed on an enhancement map
%   R = rodGranulometryEnhance(I);
%   [BW, L] = watershedSegment(R, 'method', 'marker', ...
%                                  'threshold', 0.4, 'hMinima', 0.02);
%   figure;
%   subplot(1,3,1); imshow(I,   []);            title('Raw');
%   subplot(1,3,2); imshow(BW,  []);            title('Binary mask');
%   subplot(1,3,3); imshow(label2rgb(L,'hsv','k')); title('Labels');
%
%   % Distance-transform watershed (compact objects)
%   [BW, L] = watershedSegment(R, 'method', 'distance', 'threshold', 0.5);
%
%   % Gradient watershed (sharp-boundary objects)
%   [BW, L] = watershedSegment(R, 'method', 'gradient', ...
%                                  'smoothSigma', 1.5, 'hMinima', 0.03);
%
% See also: localThresholdFast, bwdist, watershed, imextendedmax,
%           imimposemin, imhmin, label2rgb

%% ---------------- input parsing ----------------
p = inputParser;
addRequired(p, 'I', @(x) isnumeric(x) && ismatrix(x));
addParameter(p, 'method',      'marker');
addParameter(p, 'threshold',   0.5);
addParameter(p, 'threshLow',   []);
addParameter(p, 'smoothSigma', 1.0);
addParameter(p, 'hMinima',     []);    % set per-method below if empty
addParameter(p, 'minArea',     10);
addParameter(p, 'seedErosion', 2);
addParameter(p, 'bgDilation',  3);
parse(p, I, varargin{:});

method      = lower(string(p.Results.method));
thresh      = p.Results.threshold;
threshLow   = p.Results.threshLow;
sigma       = p.Results.smoothSigma;
hMin        = p.Results.hMinima;
minArea     = p.Results.minArea;
seedErosion = p.Results.seedErosion;
bgDilation  = p.Results.bgDilation;

I = im2single(I);

%% ============================================================
%% distance-transform watershed
%% ============================================================
if method == "distance"
    if isempty(hMin), hMin = 1.0; end

    BW0  = applyThreshold(I, thresh, threshLow);
    BW0  = imopen(BW0, strel('disk', 1));   % remove isolated noise pixels

    D    = bwdist(~BW0);
    Dsup = imhmin(-D, hMin);                % suppress shallow ridges
    Lraw = watershed(Dsup);
    Lraw(~BW0) = 0;                         % mask to foreground

    BW0  = bwareaopen(Lraw > 0, minArea);
    Lraw(~BW0) = 0;
    L = uint16(labelmatrix(bwconncomp(BW0)));
    BW = single(BW0);
    return
end

%% ============================================================
%% gradient-magnitude watershed
%% ============================================================
if method == "gradient"
    if isempty(hMin), hMin = 0.02; end

    Ism  = imgaussfilt(I, sigma);
    G    = imgradient(Ism);                 % gradient magnitude
    Gsup = imhmin(G, hMin);

    Lraw = watershed(Gsup);

    % Mask to foreground if threshold is provided (> 0)
    BW0  = applyThreshold(I, thresh, threshLow);
    Lraw(~BW0) = 0;
    Lraw(Lraw == 0 & BW0) = 0;            % watershed ridges inside fg → bg

    BW0  = bwareaopen(Lraw > 0, minArea);
    Lraw(~BW0) = 0;
    L  = uint16(labelmatrix(bwconncomp(BW0)));
    BW = single(BW0);
    return
end

%% ============================================================
%% marker-controlled watershed (default)
%% ============================================================
if method == "marker"
    if isempty(hMin), hMin = 0.02; end

    Ism = imgaussfilt(I, sigma);
    BW0 = applyThreshold(I, thresh, threshLow);
    BW0 = imopen(BW0, strel('disk', 1));   % clean up binary mask

    % --- foreground seeds: extended maxima restricted to foreground ------
    fgMarkers = imextendedmax(Ism, hMin) & BW0;

    % Fallback: if no extended maxima survive (very flat image), fall back
    % to eroded foreground mask as seeds
    if ~any(fgMarkers(:))
        se = strel('disk', max(1, seedErosion));
        fgMarkers = imerode(BW0, se);
    end

    % --- background marker: region safely outside foreground -------------
    bgMarker = ~imdilate(BW0, strel('disk', bgDilation));

    % --- relief image: negate so object centres are deep valleys ---------
    relief = imimposemin(-Ism, fgMarkers | bgMarker);

    % --- watershed + mask ------------------------------------------------
    Lraw = watershed(relief);
    Lraw(~BW0) = 0;

    BW0  = bwareaopen(Lraw > 0, minArea);
    Lraw(~BW0) = 0;
    L  = uint16(labelmatrix(bwconncomp(BW0)));
    BW = single(BW0);
    return
end

error('watershedSegment:unknownMethod', ...
      'Unknown method "%s". Choose ''distance'', ''gradient'', or ''marker''.', ...
      method);

end


function BW0 = applyThreshold(I, thresh, threshLow)
% applyThreshold  Binarise I with optional hysteresis.
%   Single-threshold:  BW0 = I > thresh  (when threshLow is empty).
%   Hysteresis:        seed pixels with I > thresh; grow into connected
%                      weak pixels (I > threshLow) by connected-component
%                      propagation.
if isempty(threshLow)
    BW0 = I > thresh;
else
    BW_strong = I > thresh;
    BW_weak   = I > threshLow;
    % Label weak-threshold connected components
    Lw = uint16(bwlabeln(BW_weak));
    % Keep components that contain at least one strong-threshold pixel
    keepIds = unique(Lw(BW_strong));
    keepIds(keepIds == 0) = [];
    if isempty(keepIds)
        BW0 = false(size(I));
    else
        BW0 = ismember(Lw, keepIds);
    end
end
end
