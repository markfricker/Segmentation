#!/usr/bin/env python
"""cellposeServer.py  —  persistent Cellpose worker for MATLAB integration.

Start once from PowerShell or MATLAB startup; handles all cellpose requests
for the session without per-call subprocess overhead.

Usage:
    python cellposeServer.py <watch_dir>

Watches <watch_dir> for *.req.mat files written by cellposeSegment.m.
Each request file contains the image and parameters.
Results are written as *.res.mat; errors as *.err files.

Stop the server by creating an 'exit.req' file in watch_dir, or just
close the PowerShell window / call terminate from MATLAB.
"""
import os, sys

# ---- Strip MATLAB runtime from PATH before ANY other imports ----------------
# schtasks / start /B inherit MATLAB's environment, including its runtime DLLs
# (MKL, C++ runtime, CUDA) in PATH.  cellpose.dynamics imports scipy.ndimage
# and fastremap whose C extensions find MATLAB's DLLs first and deadlock.
# Must happen before numpy/scipy are imported — before everything else.
os.environ['PATH'] = os.pathsep.join(
    p for p in os.environ.get('PATH', '').split(os.pathsep)
    if not any(s in p.upper() for s in ('MATLAB', 'MWE_', 'POLYSPACE')))
for _v in ('PYTHONPATH', 'PYTHONHOME', 'PYTHONSTARTUP'):
    os.environ.pop(_v, None)

# Numba (required by cellpose.dynamics) probes for CUDA via LLVM at import
# time.  In non-interactive Windows sessions (Task Scheduler, start /B) this
# CUDA probe hangs indefinitely.  Disabling the CUDA backend lets the CPU JIT
# initialise normally — cellpose's flow algorithms remain JIT-compiled on CPU.
os.environ['NUMBA_DISABLE_CUDA'] = '1'
# Also disable JIT compilation at import time — cellpose.dynamics triggers
# numba JIT when the module is loaded.  CPU JIT hangs in non-interactive
# sessions during LLVM IR compilation.  With JIT disabled cellpose falls back
# to interpreted Python for flow functions (slower but correct).
os.environ['NUMBA_DISABLE_JIT'] = '1'

import time, traceback
from pathlib import Path

if len(sys.argv) < 2:
    print('Usage: cellposeServer.py <watch_dir>', file=sys.stderr)
    sys.exit(1)

watch_dir = Path(sys.argv[1])
watch_dir.mkdir(parents=True, exist_ok=True)

# Write PID file so MATLAB can check if server is running
pid_file = watch_dir / 'server.pid'
pid_file.write_text(str(os.getpid()))
print(f'[server] PID {os.getpid()}  watching {watch_dir}', flush=True)

# ---- startup diagnostic log -------------------------------------------------
# Written to file because Task Scheduler output is not visible.
# Check cellpose_work/server_startup.log to see which import hangs.
_log_path = watch_dir / 'server_startup.log'
def _log(msg):
    with open(_log_path, 'a') as _f:
        import time as _t
        _f.write(f'{_t.time():.3f}  {msg}\n')
        _f.flush()

_log(f'started  PID={os.getpid()}  Python={sys.version.split()[0]}')
_log(f'PATH={os.environ.get("PATH","")[:300]}')

# ---- imports ----------------------------------------------------------------
_log('importing numpy...')
import numpy as np
_log('importing scipy.io...')
import scipy.io as sio
_log('basic imports done')

# ---- model cache ------------------------------------------------------------
_models = {}

def get_model(model_name):
    """Load model on first use, reuse thereafter."""
    if model_name not in _models:
        model_dir = Path.home() / '.cellpose' / 'models'
        mdl = model_name   # fallback: let cellpose resolve/download
        for suffix in ('', '.npy', '_0', '_1'):
            candidate = model_dir / (model_name + suffix)
            if candidate.exists():
                mdl = str(candidate)
                break
        print(f'[server] loading {mdl} (cpu) ...', flush=True)
        _log('importing torch...')
        import torch
        _log(f'torch ok: {torch.__version__}')
        _log('importing cv2...')
        import cv2; _log(f'cv2 ok: {cv2.__version__}')
        _log('importing cellpose.transforms...')
        from cellpose import transforms; _log('transforms ok')
        _log('importing scipy.ndimage...')
        import scipy.ndimage; _log('scipy.ndimage ok')
        _log('importing fastremap...')
        import fastremap; _log('fastremap ok')
        _log('importing numba (if present)...')
        try:
            import numba; _log(f'numba ok: {numba.__version__}')
        except ImportError:
            _log('numba not installed')
        _log('importing tqdm...')
        import tqdm; _log(f'tqdm ok: {tqdm.__version__}')
        _log('importing cellpose.dynamics...')
        from cellpose import dynamics; _log('dynamics ok')
        _log('importing cellpose.models.CellposeModel...')
        t0 = time.time()
        from cellpose.models import CellposeModel
        _log(f'CellposeModel ok: {time.time()-t0:.1f}s')
        # Force CPU: CUDA initialisation hangs in non-interactive Windows sessions
        # (Task Scheduler, start /B, start /MIN, etc).  GPU mode requires the
        # server to be started from an interactive PowerShell session by the user.
        # CPU performance is adequate: ~0.6s for 64px, ~6s for 256px.
        _models[model_name] = CellposeModel(pretrained_model=mdl, gpu=False)
        print(f'[server] model ready in {time.time()-t0:.1f}s', flush=True)
    return _models[model_name]

# Pre-load default model so first request is fast
try:
    get_model('cyto3')
    print('[server] ready', flush=True)
except Exception as e:
    print(f'[server] WARNING: could not pre-load cyto3: {e}', flush=True)
    print('[server] ready (model will load on first request)', flush=True)

# ---- request loop -----------------------------------------------------------
while True:
    # Exit signal
    if (watch_dir / 'exit.req').exists():
        try:
            (watch_dir / 'exit.req').unlink()
        except Exception:
            pass
        print('[server] exit requested — shutting down', flush=True)
        break

    # Process pending requests (oldest first)
    for req_file in sorted(watch_dir.glob('*.req.mat')):
        stem    = req_file.stem          # e.g. 'abc123.req'
        base    = stem.replace('.req', '')
        res_file = watch_dir / f'{base}.res.mat'
        err_file = watch_dir / f'{base}.err'

        try:
            # Retry loadmat briefly — guards against seeing the file mid-write
            # (primary guard is the .tmp→.req.mat rename in the MATLAB client).
            for _attempt in range(5):
                try:
                    mat = sio.loadmat(str(req_file))
                    break
                except Exception:
                    if _attempt == 4:
                        raise
                    time.sleep(0.1)

            def s(key, default):
                """Extract a scalar from a scipy.io.loadmat value."""
                if key not in mat:
                    return default
                v = np.squeeze(mat[key])
                return v.item() if v.ndim == 0 else default

            I       = np.squeeze(mat['I']).astype(np.float32)
            model   = str(mat['model'].flat[0]) if 'model' in mat else 'cyto3'
            diam    = s('diameter',      10.0)
            cp      = s('cellprob',       0.0)
            ft      = s('flowthreshold',  0.8)
            niter   = int(s('niter',      0))
            minsize = int(s('minsize',    0))
            timeout = s('timeout',      300.0)

            m = get_model(model)

            import threading, io, contextlib
            _res = [None, None]
            def _run():
                try:
                    kwargs = dict(diameter=diam, cellprob_threshold=cp,
                                  flow_threshold=ft, do_3D=False)
                    if niter > 0:
                        kwargs['niter'] = niter
                    # Suppress Cellpose's own stdout/stderr chatter
                    with contextlib.redirect_stdout(io.StringIO()), \
                         contextlib.redirect_stderr(io.StringIO()):
                        masks, _, _ = m.eval([I], **kwargs)
                    _res[0] = masks[0].astype(np.uint16)
                except Exception as e:
                    _res[1] = traceback.format_exc()

            t0 = time.time()
            t = threading.Thread(target=_run, daemon=True)
            t.start()
            t.join(timeout=timeout)

            if t.is_alive():
                raise RuntimeError(f'eval timed out after {timeout:.0f} s')
            if _res[1] is not None:
                raise RuntimeError(_res[1])

            L = _res[0]

            # minSize filter
            if minsize > 0 and L.max() > 0:
                from skimage import measure
                props = measure.regionprops(L)
                Lnew  = np.zeros_like(L)
                k = 0
                for prop in props:
                    if prop.area >= minsize:
                        k += 1
                        Lnew[L == prop.label] = k
                L = Lnew

            sio.savemat(str(res_file), {'L': L}, format='5')
            print(f'[server] done {time.time()-t0:.1f}s  n={int(L.max())}', flush=True)

        except Exception:
            err_file.write_text(traceback.format_exc())
            print(f'[server] ERROR processing {req_file.name}:\n{traceback.format_exc()}',
                  flush=True)
        finally:
            try:
                req_file.unlink()
            except Exception:
                pass

    time.sleep(0.05)   # 50 ms poll interval

pid_file.unlink(missing_ok=True)
