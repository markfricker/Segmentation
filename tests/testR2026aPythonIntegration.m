classdef testR2026aPythonIntegration < matlab.unittest.TestCase
%TESTR2026APYTHONINTEGRATION  Probe R2026a Python integration capabilities.
%
% PURPOSE
%   Determine which Python calling mechanism works reliably for torch /
%   cellpose / future pretrained networks on this machine.  Results inform
%   whether we can retire the server workaround and use pyrun directly.
%
% USAGE
%   cd('C:\Users\dops0035\Documents\Research\Matlab Projects\Segmentation_sandbox')
%   addpath('tests')
%
%   % Safe tests only (never hang):
%   results = runtests('testR2026aPythonIntegration', 'ExcludingTag', 'MightHang');
%   table(results)
%
%   % All tests including potentially hanging ones (run individually):
%   runtests('testR2026aPythonIntegration/testInProcess_TorchImport')
%
% TEST CATEGORIES
%   Basic           — pyenv, simple pyrun, type conversions (always safe)
%   R2026aFeatures  — new External Languages / venv API in R2026a
%   OutOfProcess    — pyrun OutOfProcess mode (should be safe with firewall rule)
%   MightHang       — InProcess pyrun with torch/cellpose (tag to skip)
%   Server          — existing server mechanism (known working baseline)
%
% INTERPRETING RESULTS
%   If MightHang tests PASS  → retire the server; use pyrun InProcess directly.
%   If OutOfProcess tests PASS but MightHang fails → use OutOfProcess pyrun.
%   If both hang             → keep server architecture.

    properties (Constant)
        PollTimeout = 20    % seconds for server round-trip test
    end

    % =====================================================================
    % Setup
    % =====================================================================
    methods (TestClassSetup)
        function reportEnvironment(tc)
            pe = pyenv();
            tc.log(matlab.unittest.Verbosity.Verbose, sprintf( ...
                'MATLAB %s  |  Python %s  |  %s  |  Mode: %s', ...
                version('-release'), pe.Version, pe.Executable, pe.ExecutionMode));
        end
    end

    % =====================================================================
    % CATEGORY 1 — Basic pyrun (always safe, no torch)
    % =====================================================================
    methods (Test, TestTags={'Basic'})

        function testBasic_PyenvVersion(tc)
            % Python must be configured and version >= 3.9.
            pe    = pyenv();
            parts = str2double(strsplit(pe.Version, '.'));
            major = parts(1);  minor = parts(2);
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('Python version: %s  (major=%d  minor=%d)', pe.Version, major, minor));
            isValid = (major > 3) || (major == 3 && minor >= 9);
            tc.verifyTrue(isValid, ...
                sprintf('Python %s < 3.9 — R2026a requires 3.9–3.13.', pe.Version));
        end

        function testBasic_PyrunSimple(tc)
            % Basic pyrun arithmetic must work.
            result = pyrun("out = 1 + 1", "out");
            tc.verifyEqual(double(result), 2.0);
        end

        function testBasic_PyrunNumpy(tc)
            % numpy import and array creation must work via pyrun.
            sz = pyrun("import numpy as np; out = np.zeros((4,4)).shape[0]", "out");
            tc.verifyEqual(double(sz), 4.0);
        end

        function testBasic_TypeConv_DoubleArray(tc)
            % MATLAB double array must round-trip through pyrun correctly.
            A   = [1.0 2.0 3.0];
            out = pyrun("import numpy as np; out = float(np.array(x).sum())", ...
                        "out", x=A);
            tc.verifyEqual(double(out), 6.0, 'AbsTol', 1e-9);
        end

        function testBasic_TypeConv_String(tc)
            % MATLAB string must be receivable in Python as a string.
            s   = "hello";
            out = pyrun("out = msg.upper()", "out", msg=s);
            tc.verifyEqual(char(out), 'HELLO');
        end

        function testBasic_TypeConv_StringArray(tc)
            % R2026a added MATLAB string array <-> Python string list conversion.
            % This is new in R2026a; skip gracefully on older versions.
            try
                words = ["alpha", "beta", "gamma"];
                n     = pyrun("out = len(words)", "out", words=words);
                tc.verifyEqual(double(n), 3.0);
                tc.log(matlab.unittest.Verbosity.Verbose, ...
                    'String array conversion: SUPPORTED (R2026a feature confirmed)');
            catch ME
                tc.log(matlab.unittest.Verbosity.Verbose, ...
                    sprintf('String array conversion: NOT supported (%s)', ME.message));
                tc.assumeFail('String array conversion not available — skipping.');
            end
        end

    end

    % =====================================================================
    % CATEGORY 2 — R2026a-specific features
    % =====================================================================
    methods (Test, TestTags={'R2026aFeatures'})

        function testR2026a_MatlabRelease(tc)
            % Confirm we are running R2026a+.
            % version('-release') returns e.g. '2026a' — extract numeric year.
            rel  = version('-release');
            year = str2double(regexp(rel, '\d+', 'match', 'once'));
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('MATLAB release: %s  (year=%d)', rel, year));
            tc.verifyGreaterThanOrEqual(year, 2026, ...
                sprintf('Expected R2026a+, got %s', rel));
        end

        function testR2026a_VenvAPIExists(tc)
            % Check R2026a Python environment management capabilities.
            % The External Languages panel is UI-only; test what is scriptable.
            pe = pyenv();
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('pyenv: Version=%s  Mode=%s  Status=%s', ...
                    pe.Version, pe.ExecutionMode, pe.Status));

            % Check pyenv properties present in R2026a
            props = properties(pe);
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('pyenv properties: %s', strjoin(props, ', ')));

            % Check py.* namespace accessible (requires Python to be loaded)
            pyrun("out = 1", "out");   % ensure Python is loaded
            try
                v = py.sys.version;
                tc.log(matlab.unittest.Verbosity.Verbose, ...
                    sprintf('py.* namespace: accessible  (py.sys.version = %s)', char(v)));
            catch ME
                tc.log(matlab.unittest.Verbosity.Verbose, ...
                    sprintf('py.* namespace: NOT accessible (%s)', ME.message));
            end

            % Check if pip programmatic install is available (R2026a feature)
            hasPipInstall = ~isempty(which('pyenv')) && ...
                any(strcmp(props, 'VirtualEnvironment'));
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('VirtualEnvironment property: %d', hasPipInstall));

            tc.verifyTrue(true);   % informational — log results, always pass
        end

        function testR2026a_RequirementsTxt(tc)
            % Verify pip is accessible (prerequisite for requirements.txt workflow).
            [s, o] = system(sprintf('"%s" -m pip --version', pyenv().Executable));
            tc.log(matlab.unittest.Verbosity.Verbose, strtrim(o));
            tc.verifyEqual(s, 0, 'pip not accessible from configured Python.');
        end

    end

    % =====================================================================
    % CATEGORY 3 — OutOfProcess mode (should be safe with firewall rule)
    % =====================================================================
    methods (Test, TestTags={'OutOfProcess'})

        function testOOP_Bridge(tc)
            % OutOfProcess bridge must respond within 10 s.
            origMode = char(pyenv().ExecutionMode);
            tc.addTeardown(@() pyenv('ExecutionMode', origMode));

            if ~strcmp(pyenv().Status, 'NotLoaded')
                terminate(pyenv());
            end
            pyenv('ExecutionMode', 'OutOfProcess');

            t0  = tic;
            out = pyrun("out = 'oop_ok'", "out");
            elapsed = toc(t0);
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('OutOfProcess hello: %.2f s', elapsed));
            tc.verifyEqual(char(out), 'oop_ok');
            tc.verifyLessThan(elapsed, 10, ...
                'OutOfProcess bridge took > 10 s — firewall may be blocking.');
        end

        function testOOP_Numpy(tc)
            % numpy must import cleanly OutOfProcess.
            if ~strcmp(pyenv().ExecutionMode, 'OutOfProcess')
                if ~strcmp(pyenv().Status, 'NotLoaded'), terminate(pyenv()); end
                pyenv('ExecutionMode', 'OutOfProcess');
            end
            out = pyrun("import numpy as np; out = np.__version__", "out");
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('numpy %s (OutOfProcess)', char(out)));
            tc.verifyNotEmpty(char(out));
        end

        function testOOP_TorchImport(tc)
            % torch import OutOfProcess — may hang if bridge is slow.
            % If this hangs: Ctrl+C, then terminate(pyenv) to recover.
            if ~strcmp(pyenv().ExecutionMode, 'OutOfProcess')
                if ~strcmp(pyenv().Status, 'NotLoaded'), terminate(pyenv()); end
                pyenv('ExecutionMode', 'OutOfProcess');
            end
            t0  = tic;
            out = pyrun("import torch; out = torch.__version__", "out");
            elapsed = toc(t0);
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('torch %s (OutOfProcess)  %.1f s', char(out), elapsed));
            tc.verifyNotEmpty(char(out));
        end

    end

    % =====================================================================
    % CATEGORY 4 — InProcess/OutOfProcess with cellpose (hang on this machine)
    %   TAG: MightHang
    %   If a test hangs: Ctrl+C, then terminate(pyenv) to recover.
    % =====================================================================
    methods (Test, TestTags={'MightHang'})

        function testInProcess_TorchImport(tc)
            % *** MAY HANG — press Ctrl+C if no response within 30 s ***
            % InProcess torch import.  If R2026a fixed the deadlock this passes.
            if ~strcmp(pyenv().ExecutionMode, 'InProcess')
                if ~strcmp(pyenv().Status, 'NotLoaded'), terminate(pyenv()); end
                pyenv('ExecutionMode', 'InProcess');
            end
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                '*** InProcess torch import — may hang; Ctrl+C to abort ***');
            t0  = tic;
            out = pyrun("import torch; out = torch.__version__", "out");
            elapsed = toc(t0);
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('torch %s InProcess  %.1f s', char(out), elapsed));
            tc.verifyNotEmpty(char(out));
        end

        function testInProcess_CUDACheck(tc)
            % *** MAY HANG ***
            % InProcess torch.cuda.is_available().
            if ~strcmp(pyenv().ExecutionMode, 'InProcess')
                if ~strcmp(pyenv().Status, 'NotLoaded'), terminate(pyenv()); end
                pyenv('ExecutionMode', 'InProcess');
            end
            t0   = tic;
            cuda = pyrun("import torch; out = torch.cuda.is_available()", "out");
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('CUDA available InProcess: %d  (%.1f s)', logical(cuda), toc(t0)));
            tc.verifyTrue(islogical(logical(cuda)) || isnumeric(double(cuda)));
        end

        function testInProcess_CellposeImport(tc)
            % *** MAY HANG — this is the critical test ***
            % If this passes, pyrun can replace the server entirely.
            if ~strcmp(pyenv().ExecutionMode, 'InProcess')
                if ~strcmp(pyenv().Status, 'NotLoaded'), terminate(pyenv()); end
                pyenv('ExecutionMode', 'InProcess');
            end
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                '*** CellposeModel import InProcess — key test; may hang ***');
            t0  = tic;
            out = pyrun([ ...
                "from cellpose.models import CellposeModel; " ...
                "out = 'imported'"], "out");
            elapsed = toc(t0);
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('CellposeModel import InProcess: %.1f s', elapsed));
            tc.verifyEqual(char(out), 'imported');
        end

        function testInProcess_CellposeInference(tc)
            % *** MAY HANG ***
            % Full cellpose inference InProcess — the gold standard test.
            if ~strcmp(pyenv().ExecutionMode, 'InProcess')
                if ~strcmp(pyenv().Status, 'NotLoaded'), terminate(pyenv()); end
                pyenv('ExecutionMode', 'InProcess');
            end
            I   = single(rand(64, 64));
            t0  = tic;
            n   = pyrun([ ...
                "from cellpose.models import CellposeModel; " ...
                "import numpy as np; " ...
                "m = CellposeModel(pretrained_model='cyto3', gpu=True); " ...
                "masks,_,_ = m.eval([np.array(img).astype(np.float32)], " ...
                "    diameter=10, cellprob_threshold=0, flow_threshold=0.8, do_3D=False); " ...
                "out = int(masks[0].max())"], ...
                "out", img=double(I));
            elapsed = toc(t0);
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('Cellpose inference InProcess: %.1f s  objects=%d', ...
                    elapsed, double(n)));
            tc.verifyClass(double(n), 'double');
        end

        function testOOP_CellposeImport(tc)
            % *** MAY HANG — numba JIT hangs OutOfProcess on this machine ***
            if ~strcmp(pyenv().ExecutionMode, 'OutOfProcess')
                if ~strcmp(pyenv().Status, 'NotLoaded'), terminate(pyenv()); end
                pyenv('ExecutionMode', 'OutOfProcess');
            end
            t0  = tic;
            pyrun("from cellpose.models import CellposeModel; out = 'ok'", "out");
            elapsed = toc(t0);
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('CellposeModel import OutOfProcess: %.1f s', elapsed));
            tc.verifyLessThan(elapsed, 60);
        end

        function testOOP_CellposeInference(tc)
            % *** MAY HANG ***
            if ~strcmp(pyenv().ExecutionMode, 'OutOfProcess')
                if ~strcmp(pyenv().Status, 'NotLoaded'), terminate(pyenv()); end
                pyenv('ExecutionMode', 'OutOfProcess');
            end
            I  = single(rand(64,64));
            t0 = tic;
            n  = pyrun([ ...
                "from cellpose.models import CellposeModel; import numpy as np; " ...
                "m = CellposeModel(pretrained_model='cyto3', gpu=True); " ...
                "masks,_,_ = m.eval([np.array(img).astype(np.float32)], " ...
                "    diameter=10, cellprob_threshold=0, flow_threshold=0.8, do_3D=False); " ...
                "out = int(masks[0].max())"], "out", img=double(I));
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('Cellpose inference OutOfProcess: %.1f s', toc(t0)));
            tc.verifyClass(double(n), 'double');
        end

    end

    % =====================================================================
    % CATEGORY 5 — Server baseline (known working)
    % =====================================================================
    methods (Test, TestTags={'Server'})

        function testServer_FileIO(tc)
            % Temp .mat file round-trip (used by server communication).
            I     = single(rand(64, 64));
            tmpIn = [tempname '.mat'];
            tc.addTeardown(@() deleteIfExists(tmpIn));
            save(tmpIn, 'I', '-v6');
            tc.verifyTrue(exist(tmpIn, 'file') == 2, '.mat file not written.');
            data = load(tmpIn);
            tc.verifySize(data.I, [64 64]);
            tc.verifyClass(data.I, 'single');
        end

        function testServer_IsRunning(tc)
            % Report whether the cellpose server is currently alive.
            workDir = fullfile(tempdir, 'cellpose_work');
            pidFile = fullfile(workDir, 'server.pid');
            alive   = cpServerAlive_(pidFile);
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('Cellpose server running: %d', alive));
            if ~alive
                tc.log(matlab.unittest.Verbosity.Verbose, ...
                    ['Server not running.  Start with:\n' ...
                     '  python cellposeServer.py "' workDir '"']);
            end
            % Not a hard failure — just informational
            tc.verifyTrue(true);
        end

        function testServer_RoundTrip(tc)
            % Full request/response cycle via server (requires server running).
            workDir = fullfile(tempdir, 'cellpose_work');
            pidFile = fullfile(workDir, 'server.pid');
            tc.assumeTrue(cpServerAlive_(pidFile), ...
                'Cellpose server not running — start it first.');

            I     = single(rand(64, 64));
            reqId = sprintf('test_%d', randi(999999));
            reqFile = fullfile(workDir, [reqId '.req.mat']);
            resFile = fullfile(workDir, [reqId '.res.mat']);
            errFile = fullfile(workDir, [reqId '.err']);
            tc.addTeardown(@() deleteIfExists(reqFile));
            tc.addTeardown(@() deleteIfExists(resFile));
            tc.addTeardown(@() deleteIfExists(errFile));

            model='cyto3'; diameter=10; cellprob=0; flowthreshold=0.8; %#ok
            niter=0; minsize=0; timeout=30; %#ok
            save(reqFile, 'I','model','diameter','cellprob','flowthreshold', ...
                 'niter','minsize','timeout', '-v6');

            t0 = tic;
            while ~exist(resFile,'file') && ~exist(errFile,'file') && toc(t0) < tc.PollTimeout
                pause(0.25);
            end
            elapsed = toc(t0);

            tc.assertFalse(exist(errFile,'file') == 2, ...
                ['Server returned error: ' filereadSafe(errFile)]);
            tc.assertTrue(exist(resFile,'file') == 2, ...
                sprintf('No result after %.0f s — server may be busy.', tc.PollTimeout));

            result = load(resFile);
            tc.log(matlab.unittest.Verbosity.Verbose, ...
                sprintf('Server round-trip: %.2f s  objects=%d', elapsed, max(result.L(:))));
            tc.verifyClass(result.L, 'uint16');
            tc.verifySize(result.L, [64 64]);
        end

    end

end

% =========================================================================
% File-local helpers
% =========================================================================

function alive = cpServerAlive_(pidFile)
    alive = false;
    if ~exist(pidFile, 'file'), return; end
    try
        pid = str2double(strtrim(fileread(pidFile)));
        if isnan(pid) || pid <= 0, return; end
        [~, o] = system(sprintf('tasklist /FI "PID eq %d" /NH 2>NUL', pid));
        alive = contains(o, num2str(pid));
    catch
    end
end

function deleteIfExists(f)
    if exist(f, 'file'), delete(f); end
end

function s = filereadSafe(f)
    try, s = fileread(f); catch, s = '(unreadable)'; end
end
