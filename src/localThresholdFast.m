function [BW, T] = localThresholdFast(I, varargin)
% localThresholdFast  Fast local/adaptive thresholding
%
% USAGE
%   BW = localThresholdFast(I)
%   BW = localThresholdFast(I, Name, Value, ...)
%   [BW, T] = localThresholdFast(...)
%
% INPUTS
%   I  : 2-D grayscale image (uint8, uint16, single, double, or logical)
%
% NAME-VALUE PAIRS
%   'method'     : thresholding method (default 'sauvola')
%                  morphological (O(N) via erosion/dilation):
%                    'midgrey'    — local mid-grey threshold
%                    'bernsen'    — classical Bernsen method
%                  mean only:
%                    'bradley'    — Bradley-Roth adaptive mean
%                  mean + standard deviation:
%                    'niblack'    — Niblack
%                    'sauvola'    — Sauvola & Pietikainen (default)
%                    'wolf'       — Wolf-Jolion
%                    'nick'       — NICK
%                    'phansalkar' — Phansalkar
%   'windowsize' : odd neighbourhood size in pixels (default 25)
%   'k'          : method sensitivity parameter (default varies by method)
%   'r'          : dynamic range parameter for sauvola/wolf/phansalkar
%                  (default: 128 for integer images, 0.5 for [0,1] images)
%   'p','q'      : phansalkar shape parameters (default p=2, q=10)
%   'contrast'   : contrast threshold for bernsen (default 15)
%   'invert'     : invert binary output (default false)
%   'precision'  : 'single' (default) or 'double'
%
% OUTPUTS
%   BW : logical binary image (true = foreground)
%   T  : threshold surface (same size as I; only computed when requested)
%
% OVERVIEW
%   Global thresholding assigns a single value T to the whole image and
%   classifies pixels as foreground where I >= T.  This fails when
%   illumination is uneven or local contrast varies across the field.
%
%   Local (adaptive) thresholding instead computes an independent threshold
%   T(x,y) for each pixel based on the statistics of a rectangular
%   neighbourhood of size windowsize x windowsize centred at (x,y).
%
%   Statistical methods (Niblack, Sauvola, etc.) express T as a function
%   of the local mean mu and standard deviation sigma:
%
%       Niblack:    T = mu + k * sigma
%       Sauvola:    T = mu * (1 + k * (sigma/R - 1))
%       Wolf:       T = mu + k * ((sigma/R - 1) * (mu - min_global))
%       NICK:       T = mu + k * sqrt(sigma^2 + mu^2)
%       Phansalkar: T = mu * (1 + p*exp(-q*mu) + k*(sigma/R - 1))
%       Bradley:    T = mu * (1 - k)
%
%   Local mu and sigma^2 are computed efficiently in O(N) using
%   separable box-filter convolutions (equivalent to integral images).
%
%   Morphological methods avoid statistics entirely.  midgrey computes
%   T(x,y) = 0.5*(local_min + local_max) using erosion and dilation with
%   a flat disk structuring element.  Bernsen extends midgrey by treating
%   low-contrast regions separately.  Both run in O(N) regardless of
%   window size when using the flat disk approximation.
%
% NOTES
%   - T is only computed when the second output argument is requested,
%     avoiding the allocation cost in production use.
%   - The r parameter default adapts to image range: 128 for integer-valued
%     images (uint8/uint16 range) and 0.5 for normalised [0,1] inputs.
%   - For fluorescence microscopy images normalised to [0,1], Sauvola with
%     k=0.34, r=0.5, windowsize=25 is a reasonable starting point.
%
% REFERENCES
%   Niblack, W. (1986). An introduction to digital image processing.
%   Prentice-Hall.
%
%   Sauvola, J. and Pietikainen, M. (2000). Adaptive document image
%   binarization. Pattern Recognition, 33(2), 225-236.
%
%   Wolf, C. and Jolion, J.-M. (2003). Extraction and recognition of
%   artificial text in multimedia documents. Pattern Analysis and
%   Applications, 6(4), 309-326.
%
%   Khurshid, K., Siddiqi, I., Faure, C. and Vincent, N. (2009).
%   Comparison of Niblack inspired binarization methods for ancient
%   documents. Proceedings of ICDAR.
%
%   Phansalkar, N., More, S., Sabale, A. and Joshi, M. (2011).
%   Adaptive local thresholding for detection of nuclei in diversity
%   stained cytology images. ICCSP.
%
%   Bradley, D. and Roth, G. (2007). Adaptive thresholding using the
%   integral image. Journal of Graphics Tools, 12(2), 13-21.
%
%   Bernsen, J. (1986). Dynamic thresholding of grey-level images.
%   Proceedings of ICPR.
%
% EXAMPLE
%   % Sauvola thresholding on a normalised fluorescence image
%   BW = localThresholdFast(I, 'method', 'sauvola', ...
%                               'windowsize', 25, 'k', 0.34, 'r', 0.5);
%
%   % Morphological Bernsen — fast, no statistics
%   BW = localThresholdFast(I, 'method', 'bernsen', 'windowsize', 31);
%
%   % Retrieve threshold surface for inspection
%   [BW, T] = localThresholdFast(I, 'method', 'niblack', 'k', -0.2);
%   figure; imshow(T, []); title('Threshold surface');
%
% See also: watershedSegment, imbinarize

%% ---------------- input parsing ----------------
p = inputParser;
addRequired(p,'I', @(x) isnumeric(x) && ismatrix(x));
addParameter(p,'method','sauvola');
addParameter(p,'windowsize',25);
addParameter(p,'k',[]);
addParameter(p,'r',[]);
addParameter(p,'p',2);
addParameter(p,'q',10);
addParameter(p,'contrast',15);
addParameter(p,'invert',false);
addParameter(p,'precision','single');
parse(p,I,varargin{:});

method = lower(string(p.Results.method));
win    = p.Results.windowsize;
k      = p.Results.k;
R      = p.Results.r;
p_par  = p.Results.p;
q_par  = p.Results.q;
C      = p.Results.contrast;
invert = logical(p.Results.invert);
prec   = validatestring(p.Results.precision,{'single','double'});

%% ---------------- convert image once ----------------
Iin  = cast(I, prec);
maxI = max(Iin(:));

if isempty(R)
    if maxI > 1
        R = 128;
    else
        R = 0.5;
    end
end

%% ============================================================
%% morphological fast path
%% ============================================================

if method == "midgrey" || method == "bernsen"
    radius = floor(win/2);
    se = strel('disk', radius, 0);
    localMin = imerode(Iin, se);
    localMax = imdilate(Iin, se);
    mg = 0.5 .* (localMin + localMax);
    if method == "midgrey"
        BW = Iin >= mg;
        if nargout > 1
            T = mg;
        end
    else  % classical bernsen
        contrast   = localMax - localMin;
        grayMid    = cast(0.5, class(Iin));
        lowContrast  = contrast < C;
        highContrast = ~lowContrast;
        BW = false(size(Iin));
        BW(lowContrast)  = mg(lowContrast) >= grayMid;
        BW(highContrast) = Iin(highContrast) >= mg(highContrast);
        if nargout > 1
            T = mg;
            T(lowContrast) = grayMid;
        end
    end
    if invert, BW = ~BW; end
    return
end

%% ============================================================
%% mean-only path (bradley)
%% ============================================================

if method == "bradley"
    if isempty(k), k = 0.15; end
    pad    = floor(win/2);
    Ipad   = padarray(Iin, [pad pad], 'symmetric', 'both');
    kernel = ones(win, win, prec) / (win*win);
    localMean = conv2(Ipad, kernel, 'valid');
    Tlocal = localMean .* (1 - k);
    BW = Iin >= Tlocal;
    if invert, BW = ~BW; end
    if nargout > 1, T = Tlocal; end
    return
end

%% ============================================================
%% mean + std path
%% ============================================================

switch method
    case "niblack",     if isempty(k), k = -0.2;  end
    case "sauvola",     if isempty(k), k =  0.34; end
    case "wolf",        if isempty(k), k =  0.5;  end
    case "nick",        if isempty(k), k = -0.1;  end
    case "phansalkar",  if isempty(k), k =  0.25; end
    otherwise
        error('localThresholdFast:unknownMethod', ...
              'Unknown method "%s".', method);
end

pad    = floor(win/2);
Ipad   = padarray(Iin, [pad pad], 'symmetric', 'both');
kernel = ones(win, win, prec) / (win*win);

localMean   = conv2(Ipad, kernel, 'valid');
localMeanSq = conv2(Ipad.^2, kernel, 'valid');
localVar    = localMeanSq - localMean.^2;
localVar(localVar < 0) = 0;
localStd = sqrt(localVar);

switch method
    case "niblack"
        Tlocal = localMean + k .* localStd;

    case "sauvola"
        Tlocal = localMean .* (1 + k .* ((localStd ./ R) - 1));

    case "wolf"
        minI   = min(Iin(:));
        Tlocal = localMean + k .* (((localStd ./ R) - 1) .* (localMean - minI));

    case "nick"
        Tlocal = localMean + k .* sqrt(localStd.^2 + localMean.^2);

    case "phansalkar"
        Tlocal = localMean .* (1 + p_par .* exp(-q_par .* localMean) ...
                                 + k .* ((localStd ./ R) - 1));
end

BW = Iin >= Tlocal;

if invert, BW = ~BW; end

if nargout > 1, T = Tlocal; end

end
