#!/usr/bin/env python3
"""Compare scientific arrays and pinned ideal execution across two result folders."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('first', type=Path)
    parser.add_argument('second', type=Path)
    parser.add_argument('--output', type=Path)
    args=parser.parse_args()
    first=json.loads((args.first/'quantum-results.json').read_text())
    second=json.loads((args.second/'quantum-results.json').read_text())
    assert first['source_sha256']==second['source_sha256']
    assert first['ideal']==second['ideal']
    checks={}
    for name in ['fixture.npz','ideal-response-matrices.npz']:
        a=np.load(args.first/name);b=np.load(args.second/name)
        assert a.files==b.files
        checks[name]={key:bool(np.array_equal(a[key],b[key])) for key in a.files}
        assert all(checks[name].values())
    for a,b in zip(first['compiled_pairs'],second['compiled_pairs']):
        for key in ['all_to_all','bidirectional_line','observed_counts','ideal_readout_probabilities','line_readout_probabilities']:
            assert a[key]==b[key]
    result={'status':'passed','scope':'Independent ideal compilation and shot execution with pinned inputs, source and seeds; exact array/count comparison. Timing and noisy runs excluded.',
            'identical_ideal_result':True,'identical_native_gate_counts_and_sample_counts':True,
            'array_equality':checks,'source_sha256':first['source_sha256'],
            'first_result_sha256':hashlib.sha256((args.first/'quantum-results.json').read_bytes()).hexdigest(),
            'second_result_sha256':hashlib.sha256((args.second/'quantum-results.json').read_bytes()).hexdigest()}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
