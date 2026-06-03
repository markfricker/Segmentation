#!/usr/bin/env python
"""cellposeRun.py  —  Cellpose inference worker, launched by cellposeLaunch.py

Reads a float32 image from <input> (.mat, variable 'I'), runs Cellpose,
writes the uint16 label matrix to <output> (.mat, variable 'L').
On any error, writes a one-line message to --errfile so MATLAB can detect
failure by polling for that file.

stdout and stderr are redirected to DEVNULL by the launcher so there is no
pipe-buffer issue; all diagnostic output goes to --errfile on failure.
"""
import os, sys, time, argparse, traceback

parser = argparse.ArgumentParser()
parser.add_argument('input')
parser.add_argument('output')
parser.add_argument('--errfile',       required=True)
parser.add_argument('--model',         default='cyto3')
parser.add_argument('--diameter',      type=float, default=10)
parser.add_argument('--cellprob',      type=float, default=0)
parser.add_argument('--flowthreshold', type=float, default=0.8)
parser.add_argument('--niter',         type=int,   default=0)
parser.add_argument('--minsize',       type=int,   default=0)
parser.add_argument('--timeout',       type=float, default=300)
args = parser.parse_args()

def write_error(msg):
    try:
        with open(args.errfile, 'w') as f:
            f.write(msg)
    except Exception:
        pass

try:
    import numpy as np
    import scipy.io as sio
    from pathlib import Path
    import threading

    # ---- load image ---------------------------------------------------------
    mat = sio.loadmat(args.input)
    I   = np.squeeze(mat['I']).astype(np.float32)

    # ---- resolve model path (bypasses Cellpose 4 name redirect) -------------
    # Try the bare name first, then common extensions used by different
    # cellpose versions.  Falls back to name string (triggers download) only
    # if no local file is found.
    model_dir = Path.home() / '.cellpose' / 'models'
    mdl = args.model  # default: let cellpose handle it (may download)
    for suffix in ('', '.npy', '_0', '_1'):
        candidate = model_dir / (args.model + suffix)
        if candidate.exists():
            mdl = str(candidate)
            break

    # ---- load model ---------------------------------------------------------
    import torch
    torch.set_num_threads(1)          # prevent torch spawning extra threads
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    from cellpose.models import CellposeModel
    m = CellposeModel(pretrained_model=mdl)

    # ---- run eval in a daemon thread with timeout ---------------------------
    _result = [None, None]

    def _run_eval():
        try:
            kwargs = dict(diameter=float(args.diameter),
                          cellprob_threshold=float(args.cellprob),
                          flow_threshold=float(args.flowthreshold),
                          do_3D=False)
            if args.niter > 0:
                kwargs['niter'] = int(args.niter)
            masks, _, _ = m.eval([I], **kwargs)
            _result[0] = masks[0].astype(np.uint16)
        except Exception as e:
            _result[1] = traceback.format_exc()

    t0 = time.time()
    _t = threading.Thread(target=_run_eval, daemon=True)
    _t.start()
    _t.join(timeout=args.timeout)

    if _t.is_alive():
        write_error(f'Cellpose eval timed out after {args.timeout:.0f} s')
        sys.exit(1)

    if _result[1] is not None:
        write_error(f'Cellpose eval failed:\n{_result[1]}')
        sys.exit(1)

    L = _result[0]

    # ---- optional minSize filter --------------------------------------------
    if args.minsize > 0 and L.max() > 0:
        from skimage import measure
        props = measure.regionprops(L)
        Lnew  = np.zeros_like(L)
        newk  = 0
        for prop in props:
            if prop.area >= args.minsize:
                newk += 1
                Lnew[L == prop.label] = newk
        L = Lnew

    # ---- save result --------------------------------------------------------
    sio.savemat(args.output, {'L': L}, format='5')

except Exception:
    write_error(traceback.format_exc())
    sys.exit(1)
