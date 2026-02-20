function [BW, T] = localThresholdFast(I, varargin)
% localThresholdFast Fast local/adaptive thresholding with separated paths
% for morphological operations
%
%   BW = localThresholdFast(I)
%   BW = localThresholdFast(I, name, value, ...)
%   [BW, T] = localThresholdFast(...)
%
%   methods:
%
%     morphological (fast, O(N)):
%       'midgrey'      local mid-grey threshold
%       'bernsen'      classical Bernsen method
%
%     mean-only:
%       'bradley'      Bradley–Roth adaptive mean
%
%     mean + standard deviation:
%       'niblack'      Niblack
%       'sauvola'      Sauvola & Pietikäinen (default)
%       'wolf'         Wolf–Jolion
%       'nick'         NICK
%       'phansalkar'   Phansalkar
%
%   input:
%     I  : 2-D grayscale image (uint8, uint16, single, double, logical)
%
%   name–value pairs:
%     'method'     : thresholding method (default 'sauvola')
%     'windowsize' : odd neighbourhood size (default 25)
%     'k'          : method parameter (varies by algorithm)
%     'r'          : dynamic range parameter (sauvola/wolf/phansalkar)
%     'p','q'      : phansalkar parameters (default p=2, q=10)
%     'contrast'   : contrast threshold for bernsen (default 15)
%     'invert'     : invert binary output (default false)
%     'precision'  : 'single' (default) or 'double'
%
%   outputs:
%     BW : logical binary image (true = foreground)
%     T  : threshold surface (only returned if requested)
%
%   performance notes:
%     - morphological methods use erosion/dilation and are O(N).
%     - statistical methods use convolution-based local statistics.
%     - threshold surface T is only computed when requested.
%
% 
%    references:
%     Niblack, W. (1986). An introduction to digital image processing.
%         Prentice-Hall.
%
%     Sauvola, J., and Pietikäinen, M. (2000). Adaptive document image
%         binarization. Pattern Recognition, 33(2), 225–236.
%
%     Wolf, C., and Jolion, J.-M. (2003). Extraction and recognition of
%         artificial text in multimedia documents. Pattern Analysis and
%         Applications, 6(4), 309–326.
%
%     Khurshid, K., Siddiqi, I., Faure, C., and Vincent, N. (2009).
%         Comparison of niblack inspired binarization methods for ancient
%         documents. Proceedings of ICDAR.
%
%     Phansalkar, N., More, S., Sabale, A., and Joshi, M. (2011).
%         Adaptive local thresholding for detection of nuclei in
%         diversity stained cytology images. ICCSP.
%
%     Bradley, D., and Roth, G. (2007). Adaptive thresholding using the
%         integral image. Journal of Graphics Tools, 12(2), 13–21.
%
%     Bernsen, J. (1986). Dynamic thresholding of grey-level images.
%         Proceedings of ICPR.-------------------------------------------------------------------------

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
Iin = cast(I, prec);
maxI = max(Iin(:));

if isempty(R)
    R = (maxI > 1) * 128 + (maxI <= 1) * 0.5;
end

%% ============================================================
%% morphological fast path
%% ============================================================

if method == "midgrey" || method == "bernsen"
    radius = floor(win/2);
    se = strel('disk', radius,0);
    localMin = imerode(Iin, se);
    localMax = imdilate(Iin, se);
    mg = 0.5 .* (localMin + localMax);
    if method == "midgrey"
        BW = Iin >= mg;
        if nargout > 1
            T = mg;
        end
    else  % classical bernsen
        contrast = localMax - localMin;
        grayMid = cast(0.5, class(Iin));
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
    if invert
        BW = ~BW;
    end
    return
end

%% ============================================================
%% mean-only path (bradley)
%% ============================================================

if method == "bradley"
    if isempty(k), k = 0.15; end
    pad = floor(win/2);
    Ipad = padarray(Iin,[pad pad],'symmetric','both');
    kernel = ones(win,win,prec) / (win*win);
    localMean = conv2(Ipad, kernel, 'valid');
    Tlocal = localMean .* (1 - k);
    BW = Iin >= Tlocal;
    if invert
        BW = ~BW;
    end
    if nargout > 1
        T = Tlocal;
    end
    return
end

%% ============================================================
%% mean + std path
%% ============================================================

switch method
    case "niblack",     if isempty(k), k = -0.2; end
    case "sauvola",     if isempty(k), k = 0.34; end
    case "wolf",        if isempty(k), k = 0.5; end
    case "nick",        if isempty(k), k = -0.1; end
    case "phansalkar",  if isempty(k), k = 0.25; end
    otherwise
        error("unknown method.");
end

pad = floor(win/2);
Ipad = padarray(Iin,[pad pad],'symmetric','both');
kernel = ones(win,win,prec) / (win*win);

localMean = conv2(Ipad, kernel, 'valid');
localMeanSq = conv2(Ipad.^2, kernel, 'valid');
localVar = localMeanSq - localMean.^2;
localVar(localVar < 0) = 0;
localStd = sqrt(localVar);

minI = min(Iin(:));

switch method
    case "niblack"
        Tlocal = localMean + k .* localStd;

    case "sauvola"
        Tlocal = localMean .* (1 + k .* ((localStd ./ R) - 1));

    case "wolf"
        Tlocal = localMean + k .* (((localStd ./ R) - 1) .* (localMean - minI));

    case "nick"
        Tlocal = localMean + k .* sqrt(localStd.^2 + localMean.^2);

    case "phansalkar"
        Tlocal = localMean .* (1 + p_par .* exp(-q_par .* localMean) ...
                                 + k .* ((localStd ./ R) - 1));
end

BW = Iin >= Tlocal;

if invert
    BW = ~BW;
end

if nargout > 1
    T = Tlocal;
end

end