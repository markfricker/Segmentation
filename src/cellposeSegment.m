function [L, BW] = cellposeSegment(I, params)
% cellposeSegment  Deep-learning cell/organelle segmenter via Cellpose
%
% USAGE
%   [L, BW] = cellposeSegment(I)
%   [L, BW] = cellposeSegment(I, params)
%
% INPUTS
%   I       - 2-D grayscale image (any numeric class); converted to
%             single [0,1] internally via im2single.
%   params  - optional struct with fields:
%            .model         = 'cyto3' % Cellpose model name or absolute path
%                                     % to a custom model file
%            .diameter      = 10      % expected object diameter in pixels
%            .cellProb      = -1      % cell probability threshold in [−6,6];
%                                     % lower values detect fainter / smaller
%                                     % objects.  0 found optimal for
%                                     % mitochondria fluorescence (sweepCellpose
%                                     % 2026-03-03); try −1 for slightly higher
%                                     % recall at cost of more false positives.
%            .flowThreshold = 0.8     % flow-error threshold in [0.1,3];
%                                     % higher values detect more objects
%                                     % at the cost of boundary precision
%            .nIter         = 0       % Cellpose gradient-descent iterations.
%                                     % 0 = use Cellpose default (~200).
%                                     % Set to 2000 for cyto3 on long thin
%                                     % structures (mitochondria, bacteria).
%                                     % Uses the Python Cellpose API directly
%                                     % (MathWorks segmentCells2D does not
%                                     % expose this parameter).
%            .normalize     = true    % API consistency; L and BW are always
%                                     % integer label / binary outputs
%            .minSize       = 0       % minimum object area in px^2.
%                                     % Objects with fewer pixels are removed
%                                     % from L and BW after segmentation.
%                                     % 0 = keep all objects (no filtering).
%                                     % ~50 recommended when comparing models
%                                     % that differ in false-positive rate.
%            .timeout       = 900     % seconds to wait for the cellpose
%                                     % server before erroring out. The
%                                     % server always runs on CPU (see
%                                     % cellposeServer.py), so full-frame
%                                     % images combined with nIter = 2000
%                                     % can take several minutes; raise this
%                                     % further for large frames.
%
% OUTPUTS
%   L       - label image, uint16, same size as I.
%             Each detected object carries a unique positive integer label.
%             Background = 0.  Use label2rgb(L) for display.
%   BW      - binary segmentation mask, single precision, same size as I.
%             Pixels belonging to detected objects are 1; background is 0.
%             Compatible with the [0,1] convention of other segmentation
%             functions such as watershedSegment and refineSegment.
%
% REQUIREMENTS
%   Medical Imaging Toolbox (MathWorks)
%   Medical Imaging Toolbox Interface for Cellpose Library Add-On
%     (install via Home > Add-Ons > search "Cellpose")
%   Cellpose Python package installed in the Python environment that
%     MATLAB's pyenv points to:
%       pip install cellpose
%
% OVERVIEW
%   Wraps the MathWorks cellpose() + segmentCells2D() API introduced in the
%   Medical Imaging Toolbox Interface for Cellpose Library add-on.  A
%   Cellpose model object is created on each call; for batch processing,
%   cache the model object externally and call segmentCells2D() directly to
%   avoid repeated model loading overhead.
%
%   MODEL SELECTION
%   'cyto3'  — Cellpose 3 super-generalist cytoplasm model (DEFAULT).
%              Trained on 9 heterogeneous datasets (Nature Methods 2025);
%              best general-purpose choice for fluorescence mito images.
%              Requires Python cellpose >= 3.x (the MATLAB add-on's
%              downloadCellposeModels does not include cyto3; the function
%              downloads it automatically via Python on first call).
%   'bact_fluor_cp3' — Cellpose 3 model trained specifically on fluorescent
%              bacterial images.  Bacteria and mitochondria share similar
%              elongated morphology, making this model a strong alternative
%              to cyto3 for organelle segmentation.  Weights are downloaded
%              automatically on first call (~80 MB).  Use nIter = 2000 for
%              best results on long thin mitochondria.
%   'cyto2'  — Cellpose 2 cytoplasm model.  Weights are in the MATLAB
%              add-on's registry (downloadCellposeModels); loaded by file
%              path, no Python download needed.  Safe fallback when
%              Python cellpose 4.x is installed (cyto3 → cpsam in 4.x).
%   'cpsam'  — Cellpose-SAM (Cellpose 4 only); SAM ViT backbone fine-tuned
%              on diverse biological cell images.  Not suitable for small
%              organelles such as mitochondria (returns blank masks with
%              ImageCellDiameter = 10–80).
%   'nuclei' — Nuclear segmentation model.
%   Custom   — Pass an absolute path to a model file trained with the
%              Cellpose GUI (human-in-the-loop mode) for best mitochondria
%              specificity after fine-tuning on ~100 annotated images.
%
%   DIAMETER PARAMETER
%   Set diameter to the expected object cross-section in pixels (not the
%   long axis).  For ~8 px wide mitochondria, diameter = 8-12 is
%   appropriate.  Cellpose internally rescales the image so that objects
%   appear ~30 px wide before inference, so diameter accuracy matters more
%   for boundary quality than for whether objects are detected at all.
%
%   SENSITIVITY PARAMETERS
%   The two main knobs for improving recall on faint or small objects:
%
%   cellProb (CellThreshold, default 0, range −6 to 6):
%     Lower values accept lower-probability detections.  Empirically,
%     0 gives the best boundary quality for mitochondria fluorescence
%     images (sweepCellpose, 2026-03-03); try −1 for slightly higher
%     recall.  Values below −3 tend to produce excessive false positives.
%
%   flowThreshold (FlowErrorThreshold, default 0.4, range 0.1 to 3):
%     Acceptable error between predicted and reconstructed flow fields.
%     Increase to 0.8–1.0 to recover fragmentary or elongated objects
%     whose flow fields are imperfect; at the cost of rougher boundaries.
%
%   OUTPUT CONVENTION
%   L is the primary output: a uint16 label matrix where each detected
%   object carries a unique positive integer, and 0 is background.
%   BW is derived from L (BW = single(L > 0)) and is provided for
%   compatibility with binary-mask pipelines.
%
% NOTES
%   - The function checks for the cellpose() function at runtime and raises
%     an informative error if the add-on is not installed.
%   - Cellpose internally normalises the input image; passing normalised
%     single [0,1] or raw uint8/uint16 both give equivalent results.
%   - For 3-D stacks use segmentCells3D() from the Medical Imaging Toolbox
%     directly; cellposeSegment processes 2-D slices only.
%   - When nIter > 0, cellposeSegment bypasses segmentCells2D entirely and
%     calls the Python CellposeModel.eval() API directly (since the MATLAB
%     wrapper does not expose niter).  This requires MATLAB R2022a+ for
%     correct numpy <-> MATLAB array conversion.  In this path the model is
%     resolved by Python directly, so cpResolveModelPath is not called.
%   - When minSize > 0, objects with fewer than minSize pixels are removed
%     from L post-hoc; labels are renumbered consecutively after filtering.
%     Applied to both code paths (MATLAB add-on and Python nIter path).
%   - The MATLAB add-on's model registry (downloadCellposeModels) only
%     covers models up to Cellpose 2 ('cyto', 'cyto2', 'nuclei', ...).
%     Newer models ('cyto3', 'bact_fluor_cp3', 'cpsam') are NOT in that
%     list and must be referenced by absolute file path.
%     cellposeSegment resolves this automatically: if cellpose() rejects a
%     model name, cpResolveModelPath uses the Python cellpose package to
%     download the weights (first call only) and returns the cached file
%     path.  When nIter > 0 the Python path handles model resolution
%     natively without cpResolveModelPath.
%
% REFERENCES
%   Stringer~C., Wang~T., Michaelos~M. \& Pachitariu~M. (2021)
%   Cellpose: a generalist algorithm for cellular segmentation.
%   Nature Methods 18:100-106.
%   https://doi.org/10.1038/s41592-020-01018-x
%     -> Original Cellpose paper and cyto/nuclei models.
%
%   Pachitariu~M. \& Stringer~C. (2022)
%   Cellpose 2.0: how to train your own model.
%   Nature Methods 19:1634-1641.
%   https://doi.org/10.1038/s41592-022-01663-4
%     -> Human-in-the-loop fine-tuning for custom models.
%
%   Stringer~C. \& Pachitariu~M. (2025)
%   Cellpose3: one-click image restoration for improved cellular
%   segmentation. Nature Methods 22:592-599.
%   https://doi.org/10.1038/s41592-025-02595-5
%     -> cyto3 model and Cellpose3 paper.
%
%   MathWorks (2025) Medical Imaging Toolbox Interface for Cellpose Library.
%   https://www.mathworks.com/help/medical-imaging/ref/cellpose.html
%     -> MATLAB wrapper used by this function.
%
% EXAMPLE
%   % Basic segmentation with default cyto3 model
%   pCP.model    = 'cyto3';
%   pCP.diameter = 10;       % ~8 px wide mitochondria
%   [L, BW] = cellposeSegment(I, pCP);
%
%   % Increase recall for faint/small mitochondria
%   pCP.cellProb      = -2;  % accept lower-probability detections
%   pCP.flowThreshold = 0.8; % tolerate imperfect flow fields
%   [L, BW] = cellposeSegment(I, pCP);
%
%   % Display results
%   figure;
%   subplot(1,3,1); imshow(I,[]);          title('Raw');
%   subplot(1,3,2); imshow(BW,[]);         title('Binary mask');
%   subplot(1,3,3); imshow(label2rgb(L,'hsv','k')); title('Labels');
%
%   % Use binary mask as weight on another enhancer
%   R_rod = rodGranulometryEnhance(I) .* BW;
%
%   % cyto3 with increased niter for elongated mitochondria (literature rec.)
%   pCP.model = 'cyto3';
%   pCP.nIter = 2000;    % uses Python API directly; not exposed by MATLAB add-on
%   [L, BW] = cellposeSegment(I, pCP);
%
%   % bact_fluor_cp3: Cellpose3 model trained on fluorescent bacteria (good for mito)
%   pCP.model = 'bact_fluor_cp3';
%   pCP.nIter = 2000;
%   [L, BW] = cellposeSegment(I, pCP);
%
%   % Custom fine-tuned model
%   pCP.model = 'C:\models\mito_finetuned';
%   [L, BW] = cellposeSegment(I, pCP);
%
%   % Refine boundaries with refineSegment
%   [~, L_refined] = refineSegment(Iraw, L, 'method','dilate','maxExpand',5);
%
% See also: watershedSegment, refineSegment, localThresholdFast

% --- defaults ---------------------------------------------------------------
% PyTorch deadlocks against MATLAB threads in InProcess mode.
% OutOfProcess runs Python in a separate process, avoiding the conflict.
% Requires the Windows Firewall to allow python.exe inbound connections.
if strcmp(pyenv().Status, 'NotLoaded')
    pyenv('ExecutionMode', 'OutOfProcess');
end

if nargin < 2, params = struct(); end
if ~isfield(params, 'model'),         params.model         = 'cyto3'; end
if ~isfield(params, 'diameter'),      params.diameter      = 6;      end
if ~isfield(params, 'cellProb'),      params.cellProb      = -1;        end
if ~isfield(params, 'flowThreshold'), params.flowThreshold = 0.8;     end
if ~isfield(params, 'nIter'),         params.nIter         = 200;       end
if ~isfield(params, 'minSize'),       params.minSize       = 0;       end
if ~isfield(params, 'normalize'),     params.normalize     = true;    end
if ~isfield(params, 'timeout'),       params.timeout       = 900;     end

% --- input validation -------------------------------------------------------
if size(I, 3) > 1
    error('cellposeSegment:badInput', ...
          'cellposeSegment: expected 2-D grayscale image, got %d-channel input.', ...
          size(I,3));
end

try
    pyrun("import cellpose");
catch pyErr
    error('cellposeSegment:notFound', ...
          ['cellposeSegment: Python cellpose package not found.\n' ...
           'Install it in your Python environment:\n' ...
           '  pip install cellpose\n' ...
           'Python error: %s'], pyErr.message);
end

% Guard against Cellpose 4, whose CellposeModel uses a SAM ViT backbone
% that is architecturally incompatible with CP3 models (cyto3, bact_fluor_cp3,
% etc.).  CP4 also silently redirects the name 'cyto3' to 'cpsam', which is
% unsuitable for small organelles.  Downgrade to keep CP3 behaviour.
cpVer = char(pyrun("import importlib.metadata; out = importlib.metadata.version('cellpose')", "out"));
cpMaj = sscanf(cpVer, '%d', 1);
if ~isempty(cpMaj) && cpMaj >= 4
    error('cellposeSegment:cp4Incompatible', [ ...
        'Cellpose %s is installed, but CP4 is not compatible with CP3\n' ...
        'models (cyto3, bact_fluor_cp3, etc.) and its default model\n' ...
        '(cpsam) returns blank masks for small organelles.\n\n' ...
        'Fix: downgrade to the last CP3 release:\n' ...
        '  pip install "cellpose<4"\n\n' ...
        'Then restart MATLAB to reload the Python environment.'], cpVer);
end

I = im2single(I);

% --- run Cellpose via persistent server --------------------------------------
% pyrun and system() both hang when PyTorch/CUDA initialises inside a process
% spawned from MATLAB.  The server runs as an independent process started
% directly from PowerShell (where Python works fine), loads the model once,
% and handles all requests via temp files.  MATLAB polls for results.
%
% START THE SERVER ONCE before using cellposeSegment:
%   In PowerShell:
%     python "...Segmentation_sandbox\src\cellposeServer.py" "C:\cellpose_work"
%   Or call cellposeServerStart() from MATLAB (uses system start /B).
[L, BW] = cpRunViaServer(I, params);
end

% ---- local helper: communicate with cellposeServer.py ----------------------
function [L, BW] = cpRunViaServer(I, params)

srcDir  = fileparts(mfilename('fullpath'));
workDir = fullfile(tempdir, 'cellpose_work');
if ~exist(workDir, 'dir'), mkdir(workDir); end

% ---- ensure server is running -----------------------------------------------
pidFile   = fullfile(workDir, 'server.pid');
readyFile = fullfile(workDir, 'server.ready');
if ~cpServerAlive(pidFile)
    % Remove a stale PID file left behind by a server that has since died,
    % so cpServerAlive sees the fresh PID written by the new server rather
    % than a dead one. Also clear any stale ready marker — the process that
    % wrote it is gone, so it must not be mistaken for the new one's status.
    if exist(pidFile, 'file'), delete(pidFile); end
    if exist(readyFile, 'file'), delete(readyFile); end
    fprintf('Starting cellpose server (first call ~15 s)...\n');
    pyExeStr = char(pyenv().Executable);   % char: exist() rejects string scalars
    if isempty(pyExeStr), pyExeStr = 'python'; end
    % Use pythonw.exe (no console window) for the server subprocess.
    % Falls back to python.exe if pythonw.exe is not found alongside it.
    [pyDir, ~, pyExt] = fileparts(pyExeStr);
    pywExe = fullfile(pyDir, ['pythonw' pyExt]);
    if ispc && exist(pywExe, 'file')
        serverPyExe = pywExe;
    else
        serverPyExe = pyExeStr;
    end
    serverScript = fullfile(srcDir, 'cellposeServer.py');
    % Windows Job Objects propagate to all child processes regardless of
    % console flags (/B, /MIN, etc).  The only reliable escape is the Task
    % Scheduler service, which runs processes in its own session outside any
    % application job object.
    taskName = 'MATLABCellposeServer';
    % Write a .ps1 helper to avoid multi-layer quoting issues with paths
    % that contain spaces (e.g. "Matlab Projects\...").  Register-ScheduledTask
    % handles spaces in -Execute/-Argument cleanly; schtasks /TR does not.
    ps1File = fullfile(workDir, 'startServer.ps1');
    fid = fopen(ps1File, 'w');
    fprintf(fid, 'Unregister-ScheduledTask -TaskName ''%s'' -Confirm:$false -ErrorAction SilentlyContinue\r\n', taskName);
    fprintf(fid, '$act = New-ScheduledTaskAction -Execute ''%s'' -Argument ''"%s" "%s"''\r\n', serverPyExe, serverScript, workDir);
    fprintf(fid, '$set = New-ScheduledTaskSettingsSet -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries -MultipleInstances IgnoreNew\r\n');
    fprintf(fid, 'Register-ScheduledTask -TaskName ''%s'' -Action $act -Settings $set -Force | Out-Null\r\n', taskName);
    fprintf(fid, 'Start-ScheduledTask -TaskName ''%s''\r\n', taskName);
    fclose(fid);
    system(sprintf('powershell -WindowStyle Hidden -NonInteractive -NoProfile -ExecutionPolicy Bypass -File "%s"', ps1File));
    % Wait for the server to actually finish loading and enter its request
    % loop (server.ready), not just for the process to start (server.pid,
    % written well before the model finishes loading) — up to 120 s.
    t0 = tic;
    while ~cpServerReady(pidFile, readyFile) && toc(t0) < 120
        pause(0.5);
    end
    if ~cpServerReady(pidFile, readyFile)
        error('cellposeSegment:serverTimeout', ...
              ['Cellpose server did not become ready within 120 s.\n' ...
               'Start it manually from PowerShell:\n' ...
               '  python "%s" "%s"'], serverScript, workDir);
    end
    fprintf('Cellpose server ready.\n');
else
    % Server process is already alive but may have JUST been started by
    % another concurrent call and still be mid-model-load — wait briefly
    % for it to finish rather than writing a request into an unattended
    % queue with no visible progress.
    t0 = tic;
    while ~cpServerReady(pidFile, readyFile) && toc(t0) < 120
        pause(0.25);
    end
end

% ---- write request ----------------------------------------------------------
reqId  = sprintf('%s_%d', datestr(now,'yyyymmddHHMMSSFFF'), randi(99999));
reqFile = fullfile(workDir, [reqId '.req.mat']);
resFile = fullfile(workDir, [reqId '.res.mat']);
errFile = fullfile(workDir, [reqId '.err']);

model         = params.model;        %#ok<NASGU>
diameter      = params.diameter;     %#ok<NASGU>
cellprob      = params.cellProb;     %#ok<NASGU>
flowthreshold = params.flowThreshold;%#ok<NASGU>
niter         = params.nIter;        %#ok<NASGU>
minsize       = params.minSize;      %#ok<NASGU>
timeout       = params.timeout;      %#ok<NASGU>
% Write to a .tmp file then rename so the server never sees a partial write.
tmpFile = [reqFile '.tmp'];
save(tmpFile, 'I', 'model', 'diameter', 'cellprob', 'flowthreshold', ...
     'niter', 'minsize', 'timeout', '-v6');
movefile(tmpFile, reqFile);

% ---- poll for result --------------------------------------------------------
% Give the client an extra 30 s over the server-side eval timeout so the
% server's own timeout error (written to errFile) has time to land first.
pollTimeout = params.timeout + 30;
t0 = tic;
while ~exist(resFile, 'file') && ~exist(errFile, 'file')
    if toc(t0) > pollTimeout
        try, delete(reqFile); catch, end
        error('cellposeSegment:timeout', 'Cellpose timed out after %d s.', pollTimeout);
    end
    pause(0.25);
end

if exist(errFile, 'file')
    msg = fileread(errFile);
    delete(errFile);
    error('cellposeSegment:failed', 'Cellpose error:\n%s', msg);
end

result = load(resFile);
delete(resFile);
L  = uint16(result.L);
BW = single(L > 0);
end

% ---- helper: is the server process alive? -----------------------------------
function alive = cpServerAlive(pidFile)
alive = false;
if ~exist(pidFile, 'file'), return; end
try
    pid = str2double(strtrim(fileread(pidFile)));
    if isnan(pid) || pid <= 0, return; end
    % Hidden PowerShell check avoids cmd window flash on every poll.
    % char(34) used for the double-quote delimiters around the -Command
    % argument — MATLAB R2017a+ parses bare "..." as a string-array literal,
    % so they cannot appear as source characters inside a char-array literal.
    dq = char(34);
    psCmd = sprintf( ...
        ['powershell -WindowStyle Hidden -NonInteractive -NoProfile -Command ' ...
         dq '(Get-Process -Id %d -ErrorAction SilentlyContinue) -ne $null' dq], pid);
    [~, out] = system(psCmd);
    alive = contains(strtrim(out), 'True');
catch
end
end

% ---- helper: is the server actually able to serve requests? -----------------
% Distinct from cpServerAlive — the PID file is written the moment the
% process starts, well before the model finishes loading (up to ~15-120 s).
% server.ready is written by cellposeServer.py only once it has entered its
% request loop, so this is the correct readiness signal to wait on before
% writing a request.
function ready = cpServerReady(pidFile, readyFile)
ready = cpServerAlive(pidFile) && exist(readyFile, 'file');
end
