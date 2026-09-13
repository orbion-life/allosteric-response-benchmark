"""Future campaign worker: fixed inherited physics with explicit action stops."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import nonlinear_profile as profile

ROOT=Path(__file__).resolve().parent


def run(protocol_path, case_index, output):
    p=json.loads(protocol_path.read_text())
    for name,expected in p['executing_source_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=expected:
            raise ValueError(f'Frozen continuation source changed: {name}')
    case=p['cases'][case_index]
    control={'query':-1,'base':0,'operator':None}
    Base=profile.CountedOperator

    class LimitedOperator(Base):
        def __init__(self,H):
            super().__init__(H);control['operator']=self

        def action(self,x,label):
            at=control['query']
            count=1 if x.ndim==1 else x.shape[1]
            prospective=self.calls['vector_equivalents']-control['base']+count
            if at>=0 and prospective>case['taylor_vector_equivalent_stops'][at]:
                raise RuntimeError(f'Taylor query{at} action budget exceeded: next={prospective}')
            return super().action(x,label)

    original=profile.expm_multiply

    def limited_exp(*args,**kwargs):
        control['query']+=1
        if control['query']>=3:
            raise RuntimeError('An unlisted Taylor query was requested.')
        control['base']=control['operator'].calls['vector_equivalents']
        return original(*args,**kwargs)

    profile.CountedOperator=LimitedOperator
    profile.expm_multiply=limited_exp
    started=time.perf_counter()
    try:
        profile.run(output,case['faces'],case['spacing'])
        raw=output/'receipt.json'
        component=json.loads(raw.read_text())
        if component['shape']!=case['shape'] or component['states']!=case['states']:
            raise AssertionError('Executed continuation grid differs from the declared case.')
        raw.rename(output/'component-receipt.json')
        receipt=dict(status='COMPLETE',continuation_protocol_sha256=hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
                     protocol_id=p['protocol_id'],case=case,
                     component_receipt='component-receipt.json',
                     component_receipt_sha256=hashlib.sha256((output/'component-receipt.json').read_bytes()).hexdigest(),
                     interpretation='The component receipt identifies the inherited physical implementation and its original protocol; this continuation protocol authorizes the new grid, timing and action caps.',
                     wrapper_elapsed_seconds=time.perf_counter()-started,
                     independent_comparison_pass=component['independent_comparison_pass'])
        (output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    except BaseException as exc:
        output.mkdir(parents=True,exist_ok=True)
        (output/'continuation-failure.json').write_text(json.dumps(dict(status='FAILED',type=type(exc).__name__,message=str(exc),
            protocol_id=p['protocol_id'],case=case,query_index=control['query'],
            operator_counts=None if control['operator'] is None else control['operator'].calls),indent=2)+'\n')
        raise


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--protocol',type=Path,required=True)
    parser.add_argument('--case-index',type=int,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();run(args.protocol,args.case_index,args.output)
