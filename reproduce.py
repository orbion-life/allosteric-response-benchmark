"""One-command reproduction, independent verification, figures and provenance."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT/"results")
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["PULSAR_OUTPUT_DIR"] = str(output)
    # Single-thread settings reduce platform-dependent scheduling variability.
    env.update({"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "MPLBACKEND": "Agg"})
    began = time.perf_counter()
    for script in ("experiment.py", "sensitivity.py", "independent_check.py"):
        print(f"Running {script}", flush=True)
        result = subprocess.run([sys.executable, "-W", "error", str(ROOT/script)], env=env, cwd=ROOT,
                                text=True, capture_output=True)
        (output/(Path(script).stem+".log")).write_text(result.stdout+result.stderr)
        if result.returncode:
            print(result.stdout+result.stderr, file=sys.stderr)
            raise SystemExit(result.returncode)
    from verify_results import verify
    receipt = verify(output)
    from make_figures import make_figures
    make_figures(output, output/"figures")
    source_files = sorted(ROOT.glob("*.py")) + [ROOT/"requirements-lock.txt"]
    provenance = {"status": "completed", "warnings_policy": "Python warnings raised as errors in all numerical subprocesses",
        "python": sys.version, "platform": platform.platform(),
        "wall_seconds": time.perf_counter()-began,
        "source_sha256": {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in source_files},
        "verification": receipt}
    (output/"execution.json").write_text(json.dumps(provenance, indent=2)+"\n")
    print(json.dumps(provenance, indent=2))

if __name__ == "__main__":
    main()
