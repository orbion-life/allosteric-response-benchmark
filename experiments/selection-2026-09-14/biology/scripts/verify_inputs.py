"""Independently restore apo coordinates and the bound-reference receiver mask."""
from pathlib import Path
import json, numpy as np
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1
P = Path(__file__).resolve().parents[1]
z = np.load(P / 'inputs/assay-domain-input.npz')
parser = MMCIFParser(QUIET=True)
apo = parser.get_structure('apo', P / 'inputs/1GFC.cif')[0]['A']
res = [apo[(' ', i, ' ')] for i in range(3, 59)]
assert ''.join(seq1(r.resname) for r in res) == ''.join(z['sequence'])
coords = np.array([r['CA'].coord for r in res])
error = float(np.max(abs(coords - z['r0'])))
assert error < 2e-6  # Bio.PDB stores coordinates as float32.
bound = parser.get_structure('bound', P / 'inputs/2VWF.cif')[0]
partner = np.array([a.coord for a in bound['B'].get_atoms() if a.element not in ['H', 'D']])
mask = []
for pos in range(1, 57):
    r = bound['A'][(' ', pos, ' ')]
    atoms = np.array([a.coord for a in r.get_atoms() if a.element not in ['H', 'D']])
    mask.append(float(np.linalg.norm(atoms[:, None] - partner[None], axis=2).min()) < 5)
assert np.array_equal(mask, z['receiver'])
print(json.dumps({'status': 'PASS', 'apo_coordinate_max_abs_error_A': error,
                  'receiver_canonical': z['canonical'][mask].tolist(),
                  'scope': 'Sequence, deposited apo coordinates and independently recomputed 2VWF/GAB2 contact mask; P212A reference and structure-quality caveats remain.'}, indent=2))
