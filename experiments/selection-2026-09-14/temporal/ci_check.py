"""Check frozen protocols and rerun the tiny pulse in an isolated source tree.

This deliberately requires no protein input, full kernel, GPU or network access.
Passing checks do not change the recorded four-protein temporal failures.
"""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import platform
import runpy
import shutil
import subprocess
import sys
import time

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parent


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main(output):
    started = time.perf_counter()
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    required = json.loads((ROOT / "COMPACT-CI.json").read_text())["required_files"]
    manifest = json.loads((ROOT / "COMPACT-PUBLICATION-MANIFEST.json").read_text())
    declared = {row["file"]: row for row in manifest["files"]}
    copied = json.loads((ROOT / "REPOSITORY-SELECTION.json").read_text())["copied_files"]
    for row in copied:
        path = ROOT / row["file"]
        assert path.stat().st_size == row["bytes"], row["file"]
        assert digest(path) == row["sha256"], row["file"]
        assert declared[row["file"]] == row, row["file"]
    for name in required:
        assert digest(ROOT / name) == declared[name]["sha256"], name
    correction = runpy.run_path(str(ROOT / "corrections/receiver-peak-indexing/correct_peaks.py"))
    correction_fixture = correction["fixture_check"]()

    def run(args, cwd, log):
        with (output / log).open("w") as handle:
            subprocess.run([sys.executable, "-W", "error", *args], cwd=cwd,
                           stdout=handle, stderr=subprocess.STDOUT, check=True)

    run([str(ROOT / "ci_protocol_check.py"), "--root", str(ROOT), "--kind", "base", "--output", str(output / "saved-protocol-tiny.json")],
        ROOT, "saved-protocol-tiny.log")
    run([str(ROOT / "ci_protocol_check.py"), "--root", str(ROOT), "--kind", "event", "--output",
         str(output / "saved-event-protocol.json")], ROOT, "saved-event-protocol.log")

    fresh = output / "fresh-source"
    fresh.mkdir()
    for name in required:
        if name.startswith("tiny-field/"):
            continue
        destination = fresh / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / name, destination)
    assert not (fresh / "inputs").exists()
    assert not list(fresh.rglob("*.npy"))
    run(["tiny_pulse.py"], fresh, "fresh-tiny-pulse.log")
    run([str(ROOT / "ci_protocol_check.py"), "--root", str(fresh), "--kind", "base", "--output", str(output / "fresh-protocol-tiny.json")],
        fresh, "fresh-protocol-tiny.log")
    run([str(ROOT / "ci_protocol_check.py"), "--root", str(fresh), "--kind", "event", "--output",
         str(output / "fresh-event-protocol.json")], fresh, "fresh-event-protocol.log")

    discrepancies = {}
    with np.load(ROOT / "tiny-field/raw.npz") as old, np.load(fresh / "tiny-field/raw.npz") as new:
        assert set(old.files) == set(new.files)
        for key in old.files:
            assert old[key].shape == new[key].shape, key
            discrepancies[key] = float(np.max(np.abs(old[key] - new[key])))
            assert discrepancies[key] <= 1e-8, (key, discrepancies[key])
    tiny = json.loads((fresh / "tiny-field/receipt.json").read_text())
    assert tiny["scientific_status"] == "PASS"
    assert tiny["max_probability_mass_error"] <= 1e-10
    assert tiny["minimum_probability"] >= -1e-12
    assert tiny["independent_forward_central_max_abs"] <= 1e-8
    receipt = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "status": "PASS",
        "scope": "Frozen source/protocol checks, archived tiny-pulse audit, fresh 64-state piecewise-field propagation, independent forward propagation, and saved/fresh tiny-array comparison. No four-protein reference-kernel or waveform replay.",
        "copied_public_file_hashes_verified": len(copied),
        "ci_protocol_adapter_sha256": digest(ROOT / "ci_protocol_check.py"),
        "grid_comparison_scope": "Only regenerated grid comparisons allow at most 32 binary64 ULPs. Historical exact-check outcomes and maximum discrepancies are retained in each protocol receipt. Saved times and kernel mappings are unchanged.",
        "fresh_tree_contains_protein_inputs": False,
        "corrected_receiver_peak_selector_fixture": correction_fixture,
        "peak_diagnostic_scope": "Original peak fields are withdrawn. CI checks the corrected selector fixture; replay of corrected protein peak metrics requires the separate waveform/correction workflow.",
        "fresh_tiny_array_max_abs_discrepancies": discrepancies,
        "fresh_tiny_receipt": tiny,
        "seconds": time.perf_counter() - started,
        "runtime": {"python": platform.python_version(), "numpy": np.__version__,
                    "scipy": scipy.__version__, "platform": platform.platform()},
        "scientific_status": "All twelve event-aligned protein pulse curves remain failures of the frozen 0.002 accuracy criterion.",
    }
    (output / "ci-receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, help="New directory outside archived evidence")
    main(parser.parse_args().output)
