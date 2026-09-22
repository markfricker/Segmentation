function [BW, level] = globalThresholdFast(I, varargin)
% globalThresholdFast  Fast histogram-based global (entropy) thresholding
%
% USAGE
%   BW = globalThresholdFast(I)
%   BW = globalThresholdFast(I, Name, Value, ...)
%   [BW, level] = globalThresholdFast(...)
%
% INPUTS
%   I  : 2-D grayscale image (uint8, uint16, single, double, or logical)
%
% NAME-VALUE PAIRS
%   'method'     : thresholding method (default 'li')
%                    'li'           — Li & Tan iterative minimum cross-entropy
%                    'kapur'        — Kapur-Sahoo-Wong maximum-entropy
%                    'triangle'     — Zack geometric triangle method
%                    'triangleOtsu' — triangle and Otsu, both computed in
%                                     log10 space on the nonzero pixels only,
%                                     level = min(triangle, otsu). Matches
%                                     the Frangi-threshold recipe used by
%                                     the Nellie organelle segmentation tool
%                                     (Lefebvre et al., Nat Methods 2025).
%   'bias'       : multiplier applied to the computed level (default 1.0).
%                  < 1 lowers the threshold (more foreground), > 1 raises it.
%                  A single global knob for pulling in dim structure without
%                  changing method.
%   'nbins'      : histogram bin count (default 256)
%   'invert'     : invert binary output (default false)
%   'precision'  : 'single' (default) or 'double' for the returned mask compare
%
% OUTPUTS
%   BW    : logical binary image (true = foreground, I > level*bias)
%   level : scalar threshold on the intensity scale of I, AFTER the bias
%           multiplier has been applied
%
% OVERVIEW
%   Both methods assign a single value to the whole image from the shape of
%   its intensity histogram; a pixel is foreground where I > level.  Unlike
%   the local methods in localThresholdFast, they do nothing for uneven
%   illumination — use them when the background is flat and the foreground
%   is a small, bright fraction of the field (the typical fluorescence
%   case), where Otsu's equal-variance assumption biases the threshold high
%   and drops faint structure.
%
%   LI (minimum cross-entropy).  Chooses the threshold t that minimises the
%   Kullback-Leibler divergence between the original image and its two-level
%   reconstruction (background pixels -> background mean, foreground pixels
%   -> foreground mean).  The 1998 paper replaces the exhaustive search of
%   the 1993 original with a fixed-point iteration: starting from the image
%   mean, repeatedly set
%
%       t_next = (mu_b - mu_f) / (log(mu_b) - log(mu_f))
%
%   where mu_b, mu_f are the mean intensities below and above the current t.
%   Converges in a handful of iterations.
%
%   KAPUR (maximum entropy).  Treats the normalised histogram below and
%   above the threshold as two probability distributions and picks the
%   threshold that maximises the sum of their Shannon entropies
%   H_b(t) + H_f(t).  Exhaustive over all bins (cheap on a 256-bin
%   histogram).  More sensitive to histogram noise than Li and inclined to
%   a low threshold when the dark background carries most of the entropy;
%   the two methods are offered together so they can be compared.
%
%   TRIANGLE (Zack et al., 1977).  Purely geometric: draws a line from the
%   histogram peak to the far end of its populated range, and picks the bin
%   with maximum perpendicular distance from that line.  Designed for the
%   sharply skewed, long-tailed histograms typical of ridge/vesselness-
%   filtered images (a tall background spike near zero with a thin bright
%   tail) -- exactly the case where Otsu's equal-variance assumption breaks
%   down and biases the threshold too high, and where Li/Kapur's entropy
%   assumptions are a poorer fit than the pure geometric argument. Falls
%   back to the opposite side of the peak if that side is empty.
%
%   TRIANGLEOTSU (log-space combination).  Ridge/vesselness responses often
%   span several orders of magnitude between the noise floor and real
%   structure, which is hard for either method to split well even after
%   'triangle' finds a reasonable answer in linear space.  Working in log10
%   space (excluding exact/near-zero pixels first, since log needs positive
%   values) spreads that range out so both triangle and Otsu see a more
%   usable shape; taking the minimum of the two is a conservative hedge --
%   whichever method wants to include more of the dim tail wins.  This is
%   the exact recipe Nellie (github.com/aelefebv/nellie) uses to threshold
%   its Frangi filter output.
%
% NOTES
%   - The histogram spans [min(I), max(I)] with nbins bins; the returned
%     level is a bin centre on that scale.  For an image normalised to
%     [0,1] the level is in [0,1].
%   - Degenerate inputs (flat image, one populated bin, a class that leaves
%     one side of the split empty) fall back to level = mean(I(:)).
%   - 0*log(0) is taken as 0 throughout.
%
% REFERENCES
%   Li, C.H. and Lee, C.K. (1993). Minimum cross entropy thresholding.
%   Pattern Recognition, 26(4), 617-625.
%
%   Li, C.H. and Tam, P.K.S. (1998). An iterative algorithm for minimum
%   cross entropy thresholding. Pattern Recognition Letters, 19(8), 771-776.
%
%   Kapur, J.N., Sahoo, P.K. and Wong, A.K.C. (1985). A new method for
%   gray-level picture thresholding using the entropy of the histogram.
%   Computer Vision, Graphics, and Image Processing, 29(3), 273-285.
%
%   Zack, G.W., Rogers, W.E. and Latt, S.A. (1977). Automatic measurement
%   of sister chromatid exchange frequency. Journal of Histochemistry and
%   Cytochemistry, 25(7), 741-753.
%
%   Lefebvre, A.E.Y.T. et al. (2025). Nellie: automated organelle
%   segmentation, tracking and hierarchical feature extraction in 2D/3D
%   live-cell microscopy. Nature Methods.
%
%   Sezgin, M. and Sankur, B. (2004). Survey over image thresholding
%   techniques and quantitative performance evaluation. Journal of
%   Electronic Imaging, 13(1), 146-165.
%
% EXAMPLE
%   % Li threshold on a normalised fluorescence image
%   BW = globalThresholdFast(I, 'method', 'li');
%
%   % Kapur, pulled in slightly to recover dim tubule
%   [BW, lvl] = globalThresholdFast(I, 'method', 'kapur', 'bias', 0.9);
%
%   % Triangle threshold on a Frangi/vesselness-enhanced ridge image
%   [BW, lvl] = globalThresholdFast(I, 'method', 'triangle');
%
%   % Nellie-style combined threshold on the same image
%   [BW, lvl] = globalThresholdFast(I, 'method', 'triangleOtsu');
%
% See also: localThresholdFast, watershedSegment, graythresh, otsuthresh

%% ---------------- input parsing ----------------
p = inputParser;
addRequired(p,'I', @(x) isnumeric(x) || islogical(x));
addParameter(p,'method','li');
addParameter(p,'bias',1.0);
addParameter(p,'nbins',256);
addParameter(p,'invert',false);
addParameter(p,'precision','single');
parse(p,I,varargin{:});

method = lower(string(p.Results.method));
bias   = double(p.Results.bias);
nbins  = max(2, round(p.Results.nbins));
invert = logical(p.Results.invert);
prec   = validatestring(p.Results.precision,{'single','double'});

%% ---------------- histogram ----------------
Iin = cast(I, prec);
v   = double(Iin(:));
v   = v(isfinite(v));

lo = min(v);
hi = max(v);
meanI = mean(v);

if isempty(v) || hi <= lo
    % flat / empty image — nothing to threshold
    level = meanI;
    if isnan(level), level = 0; end
    BW = false(size(Iin));
    if invert, BW = ~BW; end
    return
end

edges   = linspace(lo, hi, nbins + 1);
centres = (edges(1:end-1) + edges(2:end)) / 2;
counts  = histcounts(v, edges).';
centres = centres(:);

%% ---------------- method ----------------
switch method
    case "li"
        level = localLi(counts, centres, meanI, lo, hi, nbins);
    case "kapur"
        level = localKapur(counts, centres, meanI);
    case "triangle"
        level = localTriangle(counts, centres, meanI);
    case "triangleotsu"
        level = localTriangleOtsuLog(v, nbins, meanI);
    otherwise
        error('globalThresholdFast:unknownMethod', ...
              'Unknown method "%s". Use ''li'', ''kapur'', ''triangle'' or ''triangleOtsu''.', method);
end

%% ---------------- bias + binarise ----------------
level = bias * level;
level = min(max(level, lo), hi);

BW = Iin > cast(level, prec);
if invert, BW = ~BW; end

end

% =========================================================================
function t = localLi(counts, x, meanI, lo, hi, nbins)
% Li & Tam (1998) fixed-point iteration for minimum cross-entropy.
% The cross-entropy is defined for positive intensities, so shift the whole
% scale to start just above zero, solve there, and shift the result back.

    tol     = (hi - lo) / (2 * nbins);   % half a bin
    maxIter = 100;
    tiny    = eps(class(x));

    shift = 0;
    if lo <= 0
        shift = tol - lo;                 % smallest bin centre -> +tol
        x     = x + shift;
        meanI = meanI + shift;
    end

    t = meanI;                            % initial guess: image mean
    for k = 1:maxIter
        tPrev = t;
        fg = x > t;

        wB = sum(counts(~fg));
        wF = sum(counts(fg));
        if wB == 0 || wF == 0
            t = meanI;                    % split leaves one side empty
            break
        end

        muB = sum(x(~fg) .* counts(~fg)) / wB;
        muF = sum(x(fg)  .* counts(fg))  / wF;
        muB = max(muB, tiny);
        muF = max(muF, tiny);

        denom = log(muB) - log(muF);
        if abs(denom) < tiny
            break                         % means coincide — keep tPrev
        end
        t = (muB - muF) / denom;

        if ~isfinite(t)
            t = tPrev;
            break
        end
        if abs(t - tPrev) < tol
            break
        end
    end

    t = t - shift;                        % back to the original intensity scale
end

% =========================================================================
function t = localKapur(counts, x, meanI)
% Kapur, Sahoo & Wong (1985) maximum-entropy threshold.

    pmf = counts / sum(counts);
    nz  = pmf > 0;

    % Per-bin entropy contribution -p*log(p) with the 0*log0 = 0 convention
    plogp      = zeros(size(pmf));
    plogp(nz)  = -pmf(nz) .* log(pmf(nz));

    cumP      = cumsum(pmf);
    cumPlogp  = cumsum(plogp);
    totPlogp  = cumPlogp(end);

    best  = -inf;
    tIdx  = 0;
    for s = 1:numel(pmf) - 1
        Pb = cumP(s);
        Pf = 1 - Pb;
        if Pb <= 0 || Pf <= 0
            continue
        end
        Hb = log(Pb) + cumPlogp(s) / Pb;
        Hf = log(Pf) + (totPlogp - cumPlogp(s)) / Pf;
        H  = Hb + Hf;
        if H > best
            best = H;
            tIdx = s;
        end
    end

    if tIdx == 0
        t = meanI;
    else
        t = x(tIdx);
    end
end

% =========================================================================
function t = localTriangle(counts, x, meanI)
% Zack, Rogers & Latt (1977) geometric triangle threshold.
%
% Draws a line from the histogram peak to the far end of its populated
% range (whichever side of the peak is longer), then picks the bin with
% maximum perpendicular distance from that line -- the point where the
% histogram "bulges" furthest from a straight decay, which sits just past
% the tail of the background peak for the long-tailed histograms typical
% of ridge/vesselness-filtered images.

    [peakVal, peakIdx] = max(counts);
    nz = find(counts > 0);
    if isempty(nz) || numel(nz) < 2
        t = meanI;
        return
    end
    firstNZ = nz(1);
    lastNZ  = nz(end);

    % Use whichever side of the peak has the longer populated run --
    % that is the side with a tail worth splitting.
    if (peakIdx - firstNZ) >= (lastNZ - peakIdx)
        idxRange = firstNZ:peakIdx;
        endIdx   = firstNZ;
    else
        idxRange = peakIdx:lastNZ;
        endIdx   = lastNZ;
    end

    x1 = peakIdx; y1 = peakVal;
    x2 = endIdx;  y2 = counts(endIdx);
    lineLen = hypot(y2 - y1, x2 - x1);

    if lineLen == 0 || numel(idxRange) < 2
        t = meanI;
        return
    end

    idxRange = idxRange(:)';
    y0 = counts(idxRange)';
    d  = abs((y2 - y1) .* idxRange - (x2 - x1) .* y0 + x2*y1 - y2*x1) / lineLen;
    [~, k] = max(d);

    t = x(idxRange(k));
end

% =========================================================================
function t = localOtsu(counts, x)
% Standard Otsu (1979) between-class-variance threshold, vectorised over a
% precomputed histogram. Returns the bin centre maximising the between-class
% variance of a binary split at that bin.

    pmf = counts / sum(counts);
    cumP    = cumsum(pmf);
    cumMean = cumsum(pmf .* x);
    globalMean = cumMean(end);

    denom = cumP .* (1 - cumP);
    sigmaB2 = (globalMean .* cumP - cumMean).^2 ./ max(denom, eps(class(x)));
    sigmaB2(end) = -inf;   % split after the last bin is degenerate (cumP=1)

    [~, idx] = max(sigmaB2);
    t = x(idx);
end

% =========================================================================
function t = localTriangleOtsuLog(v, nbins, meanI)
% Nellie's Frangi-threshold recipe: triangle and Otsu, both computed in
% log10 space on the nonzero pixels only, then take the minimum of the two
% -- see the TRIANGLEOTSU note in the function header for the rationale.
% v is the full (non-log) sample vector already used for the outer
% function's linear histogram; this rebuilds its own histogram in log
% space since the preprocessing (exclude zero, log-transform) differs from
% every other method here.
%
% "Nonzero" is a hard v>0 test, but filters such as MATLAB's FrangiFilter2D
% leave floating-point noise residuals (~1e-7) across most of the nominal
% background rather than exact zeros, so v>0 barely excludes anything and
% a handful of near-machine-epsilon outlier pixels can stretch log10(v)'s
% range by several extra decades. localTriangle's "longer arm = tail"
% heuristic then mistakes that near-empty, spuriously long arm for the real
% signal tail. Building the log-histogram's range from the 0.5th/99.5th
% percentile of logV (not raw min/max) keeps a few outlier pixels from
% distorting that geometry -- histcounts silently drops the (very few)
% values outside the trimmed range, it does not error or clip them in.

    vPos = v(v > 0);
    if isempty(vPos)
        t = meanI;
        return
    end

    logV    = log10(vPos);
    logMean = mean(logV);
    loL = localQuantile(logV, 0.005);
    hiL = localQuantile(logV, 0.995);

    if hiL <= loL
        t = 10 ^ logMean;
        return
    end

    edgesL   = linspace(loL, hiL, nbins + 1);
    centresL = (edgesL(1:end-1) + edgesL(2:end)) / 2;
    centresL = centresL(:);
    countsL  = histcounts(logV, edgesL).';

    triLevel  = localTriangle(countsL, centresL, logMean);
    otsuLevel = localOtsu(countsL, centresL);

    t = 10 ^ min(triLevel, otsuLevel);
end

% =========================================================================
function q = localQuantile(x, p)
% Dependency-free quantile (avoids requiring Statistics and Machine
% Learning Toolbox's prctile/quantile for this one call). p in [0, 1].
    xs = sort(x(:));
    n  = numel(xs);
    idx = max(1, min(n, round(p * n)));
    q = xs(idx);
end
