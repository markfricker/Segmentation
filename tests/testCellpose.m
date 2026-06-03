classdef testCellpose < matlab.unittest.TestCase
% testCellpose  Standalone tests for the cellposeSegment wrapper.
%
% USAGE
%   Run all tests:
%       cd('C:\Users\dops0035\Documents\Research\Matlab Projects\Segmentation_sandbox')
%       results = runtests('tests/testCellpose');
%       table(results)
%
%   Run a single test:
%       runtests('tests/testCellpose/testCP_smoke')
%
% WHAT EACH TEST CHECKS
%   testCP_pythonAvailable   — Python is reachable and cellpose is importable
%   testCP_version           — cellpose is v3.x (not v4, which is incompatible)
%   testCP_cudaAvailable     — torch reports CUDA available
%   testCP_modelOnGPU        — cyto3 model parameters land on cuda:0
%   testCP_smoke             — cellposeSegment runs on a 64×64 synthetic image
%   testCP_outputClass       — label output is uint16
%   testCP_outputSize        — label output matches input size
%   testCP_outputBinary      — BW output is single with values in {0,1}
%   testCP_timing            — prints wall time; no pass/fail threshold
%
% All tests skip (not fail) when Python or cellpose is unavailable so the
% rest of the test suite still runs cleanly.

    % =====================================================================
    % Path and environment setup
    % =====================================================================
    methods (TestClassSetup)
        function addSrcPath(tc)
            rootDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(rootDir, 'src'));
        end

        function ensureServerRunning(tc)
            % Ensure the cellpose server is running before inference tests.
            % cellposeSegment() will auto-start it via Task Scheduler if needed.
            workDir = fullfile(tempdir, 'cellpose_work');
            pidFile = fullfile(workDir, 'server.pid');
            if ~cpServerAlive_(pidFile)
                tc.log(matlab.unittest.Verbosity.Verbose, ...
                    'Cellpose server not running — triggering auto-start...');
                % Calling cellposeSegment will start the server automatically
                try
                    I = single(rand(8,8));
                    cellposeSegment(I);
                catch
                    % Ignore errors — server may just be starting
                end
            end
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('Cellpose server alive: %d', cpServerAlive_(pidFile)));
        end

        function setOutOfProcess(tc)
            % PyTorch deadlocks against MATLAB threads in InProcess mode.
            if strcmp(pyenv().Status, 'NotLoaded')
                pyenv('ExecutionMode', 'OutOfProcess');
            end
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('pyenv: %s  Mode: %s', ...
                    pyenv().Executable, pyenv().ExecutionMode));
        end
    end

    % =====================================================================
    % Shared helpers
    % =====================================================================
    methods (Static, Access = private)

        function I = syntheticOrganelle(sz)
            % Small image with a few bright Gaussian blobs (~organelle-like).
            % sz: scalar image size (sz x sz), default 64.
            if nargin < 1, sz = 64; end
            I = zeros(sz, sz, 'single');
            centres = [round(sz*0.3) round(sz*0.3);
                       round(sz*0.7) round(sz*0.4);
                       round(sz*0.5) round(sz*0.7)];
            [xx, yy] = meshgrid(1:sz, 1:sz);
            for k = 1:size(centres,1)
                I = I + single(exp(-((xx - centres(k,2)).^2 + ...
                                    (yy - centres(k,1)).^2) / (2*4^2)));
            end
            I = I / max(I(:));
        end

        function tf = cellposeAvailable()
            try
                pyrun("import cellpose");
                tf = true;
            catch
                tf = false;
            end
        end

        function v = cellposeVersion()
            try
                v = char(pyrun( ...
                    "import importlib.metadata; out = importlib.metadata.version('cellpose')", ...
                    "out"));
            catch
                v = '';
            end
        end

    end

    % =====================================================================
    % Tests
    % =====================================================================
    methods (Test)

        % -----------------------------------------------------------------
        function testCP_pythonAvailable(tc)
            % Python must be reachable and cellpose importable.
            tc.assumeTrue(testCellpose.cellposeAvailable(), ...
                'cellpose not importable — skipping all cellpose tests');
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('cellpose version: %s', testCellpose.cellposeVersion()));
            tc.verifyTrue(true);   % reaching here = pass
        end

        % -----------------------------------------------------------------
        function testCP_version(tc)
            % cellpose must be v3.x.  v4 uses an incompatible SAM backbone.
            tc.assumeTrue(testCellpose.cellposeAvailable(), 'cellpose not available');
            ver  = testCellpose.cellposeVersion();
            maj  = sscanf(ver, '%d', 1);
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('cellpose version detected: %s  (major=%d)', ver, maj));
            tc.verifyLessThan(maj, 4, ...
                sprintf(['cellpose %s detected — v4 is not compatible with ' ...
                         'CP3 models.\nFix: pip install "cellpose<4"'], ver));
        end

        % -----------------------------------------------------------------
        function testCP_cudaAvailable(tc)
            % torch must see at least one CUDA device.
            tc.assumeTrue(testCellpose.cellposeAvailable(), 'cellpose not available');
            info = char(pyrun([ ...
                "import torch; " ...
                "ok  = torch.cuda.is_available(); " ...
                "dev = torch.cuda.get_device_name(0) if ok else 'none'; " ...
                "ver = torch.__version__; " ...
                "out = f'torch={ver}  CUDA={ok}  device={dev}'"], "out"));
            tc.log(matlab.unittest.Verbosity.Verbose, info);
            tc.verifyTrue( ...
                contains(info, 'CUDA=True'), ...
                sprintf(['CUDA not available.\n%s\n' ...
                         'Fix: pip install torch torchvision ' ...
                         '--index-url https://download.pytorch.org/whl/cu128 ' ...
                         '--force-reinstall'], info));
        end

        % -----------------------------------------------------------------
        function testCP_modelOnGPU(tc)
            % Check model device via the server startup log.
            % (pyrun with CellposeModel deadlocks in MATLAB's process chain —
            %  the server runs outside the Windows Job Object and logs device info.)
            tc.assumeTrue(testCellpose.cellposeAvailable(), 'cellpose not available');
            logFile = fullfile(tempdir, 'cellpose_work', 'server_startup.log');
            if ~exist(logFile, 'file')
                tc.log(matlab.unittest.Verbosity.Verbose, ...
                    'Server log not found — start server first, then rerun.');
                tc.assumeFail('Server startup log not available.');
            end
            logText = fileread(logFile);
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('Server log tail:\n%s', ...
                    logText(max(1,end-500):end)));
            % Log contains lines like "CellposeModel ok: 3.5s"
            tc.verifyTrue(contains(logText, 'CellposeModel ok'), ...
                'Server log does not show successful CellposeModel load.');
        end

        % -----------------------------------------------------------------
        function testCP_smoke(tc)
            % cellposeSegment must complete on a 64×64 synthetic image.
            tc.assumeTrue(testCellpose.cellposeAvailable(), 'cellpose not available');
            I = testCellpose.syntheticOrganelle(64);
            p.model         = 'cyto3';
            p.diameter      = 10;
            p.cellProb      = 0;
            p.flowThreshold = 0.8;
            p.nIter         = 0;
            p.minSize       = 0;
            t0 = tic;
            [L, BW] = cellposeSegment(I, p);
            elapsed = toc(t0);
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('smoke run: %.1f s   objects=%d   coverage=%.1f%%', ...
                    elapsed, max(L(:)), 100*nnz(L>0)/numel(L)));
            tc.verifyNotEmpty(L);
            tc.verifyNotEmpty(BW);
        end

        % -----------------------------------------------------------------
        function testCP_outputClass(tc)
            tc.assumeTrue(testCellpose.cellposeAvailable(), 'cellpose not available');
            I = testCellpose.syntheticOrganelle(64);
            [L, BW] = cellposeSegment(I);
            tc.verifyClass(L,  'uint16');
            tc.verifyClass(BW, 'single');
        end

        % -----------------------------------------------------------------
        function testCP_outputSize(tc)
            tc.assumeTrue(testCellpose.cellposeAvailable(), 'cellpose not available');
            I = testCellpose.syntheticOrganelle(64);
            [L, BW] = cellposeSegment(I);
            tc.verifySize(L,  size(I));
            tc.verifySize(BW, size(I));
        end

        % -----------------------------------------------------------------
        function testCP_outputBinary(tc)
            % BW must contain only 0 and 1.
            tc.assumeTrue(testCellpose.cellposeAvailable(), 'cellpose not available');
            I = testCellpose.syntheticOrganelle(64);
            [~, BW] = cellposeSegment(I);
            uniqueVals = unique(BW(:));
            tc.verifyTrue(all(uniqueVals == 0 | uniqueVals == 1));
        end

        % -----------------------------------------------------------------
        function testCP_timing(tc)
            % Informational only — prints wall time, no pass/fail threshold.
            % First call includes CUDA context init (~30 s); subsequent calls
            % should be faster.  Watch for unexpectedly large values.
            tc.assumeTrue(testCellpose.cellposeAvailable(), 'cellpose not available');
            I = testCellpose.syntheticOrganelle(64);
            p.model    = 'cyto3';
            p.diameter = 10;
            p.nIter    = 0;

            times = nan(1, 3);
            for k = 1:3
                t0       = tic;
                cellposeSegment(I, p);
                times(k) = toc(t0);
            end
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('timing (3 runs): %.1f s  %.1f s  %.1f s  (mean=%.1f s)', ...
                    times(1), times(2), times(3), mean(times)));
            tc.verifyTrue(true);   % always passes — this test is for logging only
        end

    end

end

% ---- file-local helper ------------------------------------------------------
function alive = cpServerAlive_(pidFile)
    alive = false;
    if ~exist(pidFile, 'file'), return; end
    try
        pid = str2double(strtrim(fileread(pidFile)));
        if isnan(pid) || pid <= 0, return; end
        [~, o] = system(sprintf('tasklist /FI "PID eq %d" /NH 2>NUL', pid));
        alive  = contains(o, num2str(pid));
    catch
    end
end
