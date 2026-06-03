#!/usr/bin/env python
"""cellposeLaunch.py  — fire-and-forget launcher called by cellposeSegment.m

Starts cellposeRun.py in a completely fresh subprocess with a clean
environment (no MATLAB paths, no PYTHONPATH/PYTHONHOME overrides), then
exits immediately so MATLAB's system() returns at once.

MATLAB polls for the output .mat file rather than waiting synchronously,
which avoids all pipe-buffer and blocking-call issues.
"""
import subprocess, sys, os, argparse

parser = argparse.ArgumentParser()
parser.add_argument('input')
parser.add_argument('output')
parser.add_argument('--errfile',       required=True)
parser.add_argument('--model',         default='cyto3')
parser.add_argument('--diameter',      default='10')
parser.add_argument('--cellprob',      default='0')
parser.add_argument('--flowthreshold', default='0.8')
parser.add_argument('--niter',         default='0')
parser.add_argument('--minsize',       default='0')
parser.add_argument('--timeout',       default='300')
args = parser.parse_args()

# ---- build a clean environment ----------------------------------------------
# Strip MATLAB-injected variables that override Python's own paths.
clean = {}
for k, v in os.environ.items():
    ku = k.upper()
    # Drop variables that corrupt the Python environment
    if ku in ('PYTHONPATH', 'PYTHONHOME', 'PYTHONSTARTUP', 'PYTHONUSERBASE'):
        continue
    clean[k] = v

# Remove MATLAB's bin directories from PATH to prevent DLL conflicts
path_parts = clean.get('PATH', '').split(os.pathsep)
clean['PATH'] = os.pathsep.join(
    p for p in path_parts
    if not any(s in p.upper() for s in ('MATLAB', 'MWE_', 'POLYSPACE')))

# ---- locate cellposeRun.py --------------------------------------------------
script_dir   = os.path.dirname(os.path.abspath(__file__))
run_script   = os.path.join(script_dir, 'cellposeRun.py')

# ---- launch cellposeRun.py — no wait ----------------------------------------
cmd = [sys.executable, run_script,
       args.input, args.output,
       '--errfile',       args.errfile,
       '--model',         args.model,
       '--diameter',      args.diameter,
       '--cellprob',      args.cellprob,
       '--flowthreshold', args.flowthreshold,
       '--niter',         args.niter,
       '--minsize',       args.minsize,
       '--timeout',       args.timeout]

flags = getattr(subprocess, 'CREATE_NO_WINDOW', 0)   # hide console on Windows
subprocess.Popen(cmd, env=clean, stdout=subprocess.DEVNULL,
                 stderr=subprocess.DEVNULL, creationflags=flags)

# Exit immediately — MATLAB polls for output/error files
print('launched', flush=True)
