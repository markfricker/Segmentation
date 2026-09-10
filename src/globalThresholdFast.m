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
%                    'li'    — Li & Tan iterative minimum cross-entropy
%                    'kapur' — Kapur-Sahoo-Wong maximum-entropy
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
    otherwise
        error('globalThresholdFast:unknownMethod', ...
              'Unknown method "%s". Use ''li'' or ''kapur''.', method);
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
