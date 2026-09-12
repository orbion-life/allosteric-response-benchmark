"""Show that a linear internal chart does not globally remove finite rotations.

This is an algebraic geometry check on one triangle. It computes no response,
protein score, physical population, or trajectory.
"""
from pathlib import Path
import hashlib
import json

import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parent


def main():
    points = np.array([[0.,0.,0.],[4.,0.,0.],[1.,3.5,0.]])
    edges = [(0,1),(0,2),(1,2)]
    H = np.zeros((9,9))
    for i,j in edges:
        a = points[i]-points[j]
        T = np.zeros((3,9))
        T[:,3*i:3*i+3] = np.eye(3)
        T[:,3*j:3*j+3] = -np.eye(3)
        H += T.T @ np.outer(a,a) @ T / np.dot(a,a)
    values,B = eigh(H)
    B = B[:,values>1e-8]
    values = values[values>1e-8]
    assert B.shape == (9,3)
    q = B.T @ (-2*(points-points.mean(axis=0)).ravel())
    rotated = points + (B@q).reshape(3,3)
    U = sum((np.sum((rotated[i]-rotated[j])**2)-np.sum((points[i]-points[j])**2))**2 /
            (8*np.sum((points[i]-points[j])**2)) for i,j in edges)
    error = float(np.max(np.abs(rotated-(2*points.mean(axis=0)-points))))
    assert error < 1e-12 and U < 1e-20
    result = {"status":"algebraic coordinate-domain diagnostic; no response calculation",
              "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "geometry_A":points.tolist(),"full_internal_dimensions":3,
              "positive_harmonic_eigenvalues":values.tolist(),
              "q_finite_pi_rotation":q.tolist(),"displacement_norm_A":float(np.linalg.norm(q)),
              "maximum_reconstruction_error_A":error,"U_at_finite_pi_rotation":float(U),
              "native_box_half_width_at_kappa1_4sigma_A":float(4/np.sqrt(values[0])),
              "native_box_half_width_at_kappa100_4sigma_A":float(4/np.sqrt(100*values[0])),
              "interpretation":"The linear internal chart contains a distant rigidly equivalent zero-energy configuration. Local4sigma/5sigma refinement alone cannot certify integration over allR^3 or a global nonlinear quotient of rigid rotations."}
    (ROOT/'coordinate-domain-check.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({"U":float(U),"reconstruction_error_A":error,"displacement_norm_A":float(np.linalg.norm(q))}))


if __name__ == '__main__':
    main()
