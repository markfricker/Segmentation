function results = compareCellposeBBox(I, margins)
%COMPARECELLPOSEBBOX  Time + quality comparison: full-frame vs bbox Cellpose.
%
%   results = compareCellposeBBox()
%   results = compareCellposeBBox(I)
%   results = compareCellposeBBox(I, margins)
%
%   Runs the full-frame cellposeSegment (the original approach) and the
%   bounding-box accelerated cellposeSegmentBBox at one or more bbox margins
%   on the SAME image, and reports for each:
%       - wall-clock time
%       - object count
%       - foreground coverage (%)
%       - agreement with the full-frame result (Dice of the binary masks,
%         and per-object count difference)
%
%   The full-frame run is the reference ("ground truth" for agreement); a
%   bbox run that segments identically would score Dice = 1.000 and the same
%   object count.  Lower Dice / different counts quantify how much the bbox
%   pre-segmentation changes the result.
%
%   INPUTS
%     I       : 2-D grayscale image.  Default: a synthetic 400x400 sparse
%               field of elongated organelle-like blobs (~8% coverage).
%     margins : vector of bboxMargin values to test (px).  Default [10 20 40].
%
%   OUTPUT
%     results : struct array, one element per run (full + each margin), with
%               fields .name .time .nObj .coverage .diceVsFull .countDiff.
%
%   NOTE  Cellpose runs on CPU in this setup, so absolute times are dominated
%   by image / crop area.  The bbox approach is faster only when total crop
%   area << full-frame area (i.e. sparse images); on dense images it can be
%   SLOWER due to per-crop model-eval overhead.
%
%   See also: cellposeSegment, cellposeSegmentBBox

    rootDir = fullfile(fileparts(mfilename('fullpath')), '..');
    addpath(fullfile(rootDir, 'src'));

    if nargin < 1 || isempty(I), I = makeSparseField(400, 22, 7); end
    if nargin < 2 || isempty(margins), margins = [10 20 40]; end
    I = im2single(I);

    % Shared Cellpose parameters (match the panel/dispatcher defaults).
    base = struct('model','cyto3', 'diameter',10, 'cellProb',0, ...
                  'flowThreshold',0.8, 'nIter',200, 'minSize',0, 'normalize',true);

    % Warm up the server / model so the first timed run is not penalised by
    % the one-off model load.
    fprintf('Warming up Cellpose server...\n');
    cellposeSegment(single(rand(16,16)), base);

    % ---- reference: full-frame (original approach) --------------------------
    fprintf('Running full-frame cellposeSegment (%dx%d)...\n', size(I,1), size(I,2));
    t0    = tic;
    Lfull = cellposeSegment(I, base);
    tFull = toc(t0);

    BWfull = Lfull > 0;
    results(1) = makeRow('full-frame', tFull, Lfull, BWfull, BWfull);

    % ---- bbox runs at each margin ------------------------------------------
    for m = margins(:)'
        p = base;
        p.bboxMargin = m;
        fprintf('Running cellposeSegmentBBox (margin=%d)...\n', m);
        t0 = tic;
        Lb = cellposeSegmentBBox(I, p);
        tb = toc(t0);
        results(end+1) = makeRow(sprintf('bbox m=%d', m), tb, Lb, Lb>0, BWfull); %#ok<AGROW>
    end

    % ---- report -------------------------------------------------------------
    printReport(results, tFull);
end

% =========================================================================
% helpers
% =========================================================================

function row = makeRow(name, t, L, BW, BWref)
    row.name       = name;
    row.time       = t;
    row.nObj       = double(max(L(:)));
    row.coverage   = 100 * nnz(BW) / numel(BW);
    inter          = nnz(BW & BWref);
    row.diceVsFull = 2*inter / (nnz(BW) + nnz(BWref) + eps);
end

function printReport(results, tFull)
    fprintf('\n================ Cellpose: full-frame vs bbox ================\n');
    fprintf('%-12s %8s %7s %9s %9s %8s\n', ...
        'run', 'time(s)', 'nObj', 'cov(%)', 'Dice', 'speedup');
    for k = 1:numel(results)
        r = results(k);
        if k == 1
            spd = '  ref';
        else
            spd = sprintf('%5.2fx', tFull / max(r.time, eps));
        end
        fprintf('%-12s %8.2f %7d %9.2f %9.3f %8s\n', ...
            r.name, r.time, r.nObj, r.coverage, r.diceVsFull, spd);
    end
    fprintf('==============================================================\n');
    fprintf(['Dice = overlap of binary masks vs full-frame (1.000 = identical).\n' ...
             'speedup = full-frame time / this run time (>1 = bbox faster).\n']);
end

function I = makeSparseField(sz, nBlob, halfLen)
%MAKESPARSEFIELD  Synthetic sparse field of elongated organelle-like blobs.
%   sz      : image size (sz x sz)
%   nBlob   : number of blobs
%   halfLen : nominal blob half-length (px)
    rng(7);                                   % reproducible layout
    I = zeros(sz, sz, 'single');
    [xx, yy] = meshgrid(1:sz, 1:sz);
    margin = round(0.08 * sz);
    for k = 1:nBlob
        cx  = randi([margin sz-margin]);
        cy  = randi([margin sz-margin]);
        ang = rand * pi;
        a   = halfLen * (0.7 + 0.6*rand);     % long axis
        b   = 2.2;                            % short axis (~thin tubule)
        xr  =  (xx-cx)*cos(ang) + (yy-cy)*sin(ang);
        yr  = -(xx-cx)*sin(ang) + (yy-cy)*cos(ang);
        I   = I + single(exp(-(xr.^2/(2*a^2) + yr.^2/(2*b^2))));
    end
    I = I / max(I(:));
end
