"""Independent scientific assertions on a completed one-command run."""
import json
from pathlib import Path
import numpy as np
import pytest
from verify_results import verify

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/"results"

@pytest.fixture(scope="module")
def data():
    assert (OUTPUT/"sensitivity-results.json").exists(), "Run python reproduce.py before pytest"
    return json.loads((OUTPUT/"sensitivity-results.json").read_text())

def test_archived_scientific_replay():
    assert verify(OUTPUT)["status"]=="passed"

def test_response_symmetry_and_negative_semidefiniteness(data):
    # These follow from this reversible step-response model, not biological directionality.
    for case in data["strength"]+data["grid_kappa100"]:
        r=np.array(case["R"])
        assert np.max(np.abs(r-r.T))<1e-10
        assert np.linalg.eigvalsh((r+r.T)/2).max()<1e-9

def test_independent_field_derivative():
    d=json.loads((OUTPUT/"independent-response-check.json").read_text())
    assert len(d["cases"])==6
    assert max(c["frechet_C_max_error"] for c in d["cases"])<1e-10

def test_reported_model_difference_and_refinement(data):
    def c(model,n):
        return next(v["C43"] for v in data["grid_kappa100"] if v["model"]==model and v["n"]==n)
    assert c("biquadratic",32)==pytest.approx(-.4583321803090694,abs=1e-8)
    assert c("harmonic",32)==pytest.approx(-.785145265762724,abs=1e-8)
    assert abs(c("biquadratic",64)-c("biquadratic",32))<.003
    assert abs(c("biquadratic",8)-c("biquadratic",64))>.08

def test_figure_artifacts_exist():
    for name in ("overview", "convergence"):
        for extension in ("svg", "png"):
            assert (OUTPUT/"figures"/f"{name}.{extension}").stat().st_size>1000
