"""Compare every scientific numeric field with the archived reference.

Timings, source hashes, software/platform descriptions and resource counters are
excluded explicitly. Absolute tolerance 1e-8, relative tolerance 1e-8; these are
replay tolerances, not estimates of model, grid or compression accuracy.
"""
import json
from pathlib import Path
import numpy as np

FILES = ("results.json", "sensitivity-results.json", "independent-response-check.json")
IGNORE = {"timestamp_unix", "source_sha256", "software", "resources", "execution_provenance",
          "seconds", "total_seconds", "wall_time_seconds", "max_rss_bytes"}

def compare(actual, expected, path="root"):
    count = 0
    if isinstance(expected, dict):
        for key, value in expected.items():
            if key in IGNORE:
                continue
            if key not in actual:
                raise AssertionError(f"Missing {path}.{key}")
            count += compare(actual[key], value, f"{path}.{key}")
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise AssertionError(f"Length mismatch at {path}")
        for i, (a, e) in enumerate(zip(actual, expected)):
            count += compare(a, e, f"{path}[{i}]")
    elif isinstance(expected, (int, float)) and not isinstance(expected, bool):
        if not np.isfinite(actual) or not np.isclose(actual, expected, atol=1e-8, rtol=1e-8):
            raise AssertionError(f"Numeric mismatch at {path}: {actual} versus {expected}")
        count += 1
    # Descriptive text is documentation, not a numerical replay assertion.
    return count

def verify(output):
    reference = Path(__file__).resolve().parent / "reference"
    counts = {}
    for name in FILES:
        counts[name] = compare(json.loads((output/name).read_text()), json.loads((reference/name).read_text()))
    independent = json.loads((output/FILES[-1]).read_text())
    for case in independent["cases"]:
        assert case["frechet_C_max_error"] < 1e-10
    receipt = {"status": "passed", "absolute_tolerance": 1e-8, "relative_tolerance": 1e-8,
               "numeric_fields_checked": counts,
               "scope": "scientific replay agreement, not continuum or biological accuracy"}
    (output/"verification.json").write_text(json.dumps(receipt, indent=2)+"\n")
    return receipt
