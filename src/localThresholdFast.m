function [BW, T] = localThresholdFast(I, varargin)
%LOCALTHRESHOLDFAST Fast local/adaptive image thresholding (multiple methods)
%
%   BW = LOCALTHRESHOLDFAST(I)
%   BW = LOCALTHRESHOLDFAST(I, Name, Value, ...)
%   [BW, T] = LOCALTHRESHOLDFAST(...)
%
%   Supported methods:
%     'sauvola'     Sauvola & Pietikäinen (default)
%     'niblack'     Niblack
%     'wolf'        Wolf–Jolion
%     'nick'        NICK
%     'phansalkar'  Phansalkar
%     'bradley'     Bradley–Roth
%     'midgrey'     Local mid-grey (morphological)
%     'bernsen'     Bernsen (morphological)
%
%   INPUT
%     I : 2-D grayscale image (uint8, uint16, single, double, logical)
%
%   NAME–VALUE PAIRS
%     'Method'     : Thresholding method (default 'sauvola')
%     'WindowSize' : Odd integer neighborhood size (default 25)
%     'K'          : Method parameter (see references)
%     'R'          : Dynamic range (Sauvola/Wolf/Phansalkar)
%     'P','Q'      : Phansalkar parameters (default P=2, Q=10)
%     'Contrast'   : Bernsen contrast threshold (default 15)
%     'Invert'     : Logical flag to invert binary output
%     'Precision'  : 'double' (default) or 'single'
%
%   OUTPUTS
%     BW : Logical binary image (true = foreground)
%     T  : Threshold surface
%
%   REFERENCES
%     [1] Niblack, W. (1986). Digital Image Processing.
%     [2] Sauvola & Pietikäinen (2000). Pattern Recognition.
%     [3] Wolf & Jolion (2003). Pattern Analysis & Applications.
%     [4] Khurshid et al. (2009). ICDAR. (NICK)
%     [5] Phansalkar et al. (2011). ICCSP.
%     [6] Bradley & Roth (2007). Journal of Graphics Tools.
%     [7] Bernsen (1986). ICPR.
%
% -------------------------------------------------------------------------

%% ---------------- Input parsing ----------------
p = inputParser;
addRequired(p,'I', @(x) isnumeric(x) && ndims(x)==2);
addParameter(p,'Method','sauvola');
addParameter(p,'WindowSize',25);
addParameter(p,'K',[]);
addParameter(p,'R',[]);
addParameter(p,'P',2);
addParameter(p,'Q',10);
addParameter(p,'Contrast',15);
addParameter(p,'Invert',false);
addParameter(p,'Precision','double');
parse(p,I,varargin{:});

method = lower(string(p.Results.Method));
win    = p.Results.WindowSize;
k      = p.Results.K;
R      = p.Results.R;
p_par  = p.Results.P;
q_par  = p.Results.Q;
C      = p.Results.Contrast;
invert = logical(p.Results.Invert);
prec   = validatestring(p.Results.Precision,{'double','single'});

%% ---------------- Image conversion ----------------
if strcmp(prec,'single')
    Iin = single(I);
else
    Iin = double(I);
end

maxI = max(Iin(:));
if isempty(R)
    R = (maxI > 1) * 128 + (maxI <= 1) * 0.5;
end

%% ---------------- Default parameters ----------------
switch method
    case "niblack",     if isempty(k), k = -0.2; end
    case "sauvola",     if isempty(k), k = 0.34; end
    case "wolf",        if isempty(k), k = 0.5; end
    case "nick",        if isempty(k), k = -0.1; end
    case "phansalkar",  if isempty(k), k = 0.25; end
    case "bradley",     if isempty(k), k = 0.15; end
    case {"midgrey","bernsen"}
    otherwise
        error("Unknown method.");
end

%% ---------------- Compute local statistics ----------------
needsStd = ismember(method,["niblack","sauvola","wolf","nick","phansalkar"]);

if needsStd
    pad = floor(win/2);
    Ipad = padarray(Iin,[pad pad],'symmetric','both');
    kernel = ones(win,win,class(Ipad)) / (win*win);

    localMean = conv2(Ipad, kernel, 'valid');
    localMeanSq = conv2(Ipad.^2, kernel, 'valid');
    localVar = localMeanSq - localMean.^2;
    localVar(localVar < 0) = 0;
    localStd = sqrt(localVar);
end

%% ---------------- Morphological min/max (FAST) ----------------
needsMinMax = ismember(method,["midgrey","bernsen"]);
if needsMinMax
    radius = floor(win/2);
    se = strel('disk', radius, 0);   % use 'square' for even faster if desired
    localMin = imerode(Iin, se);
    localMax = imdilate(Iin, se);
end

%% ---------------- Threshold computation ----------------
minI = min(Iin(:));

switch method
    case "niblack"
        T = localMean + k .* localStd;

    case "sauvola"
        T = localMean .* (1 + k .* ((localStd ./ R) - 1));

    case "wolf"
        T = localMean + k .* (((localStd ./ R) - 1) .* (localMean - minI));

    case "nick"
        T = localMean + k .* sqrt(localStd.^2 + localMean.^2);

    case "phansalkar"
        T = localMean .* (1 + p_par .* exp(-q_par .* localMean) ...
                            + k .* ((localStd ./ R) - 1));

    case "bradley"
        T = localMean .* (1 - k);

    case "midgrey"
        T = 0.5 .* (localMax + localMin);

    case "bernsen"
        contrast = localMax - localMin;
        T = 0.5 .* (localMax + localMin);
        T(contrast < C) = mean(Iin(:));
end

%% ---------------- Binary output ----------------
BW = Iin >= T;
if invert
    BW = ~BW;
end

end