#!/usr/bin/env python3
"""Check a proposed cohort's declarations and file hashes, not scientific validity."""
import argparse, hashlib, json
from pathlib import Path

REQUIRED = ['family_id','relatedness_evidence','source_accession','primary_evidence_doi',
            'input_structure','input_sequence_map','input_receiver_source',
            'reference_manifest_sha256','custodian','reviewer','overlap_review',
            'assembly_state_review','atom_coverage_review','numerical_admission_receipt',
            'candidate_parity_receipt','resource_admission_receipt','exclusion_log']

def audit(document, root):
    errors=[]; rows=document.get('families',[])
    if not rows: errors.append('No eligible families inventoried; cohort cannot launch.')
    if document.get('status')!='approved_before_label_access': errors.append('Protocol approval before label access is not recorded.')
    for key in ['analysis_owner','committed_effort_record','backup_owner','precision_design','sampling_rule','primary_comparator','analysis_plan']:
        if not document.get(key): errors.append(f'Missing cohort declaration: {key}')
    ids=[]
    for i,row in enumerate(rows):
        for key in REQUIRED:
            if not row.get(key): errors.append(f'Family row {i}: missing {key}')
        ids.append(row.get('family_id'))
        if row.get('development_overlap_status')!='reviewed_and_excluded_overlap': errors.append(f'Family row {i}: overlap review incomplete.')
        if row.get('minimum_candidate_count',0)<5: errors.append(f'Family row {i}: fewer than five candidates.')
        for record in row.get('files',[]):
            p=(root/record['path']).resolve()
            if not p.is_relative_to(root.resolve()): errors.append(f'Family row {i}: file outside inventory directory.'); continue
            if not p.is_file(): errors.append(f'Family row {i}: missing file {record["path"]}.');continue
            if hashlib.sha256(p.read_bytes()).hexdigest()!=record.get('sha256'): errors.append(f'Family row {i}: hash mismatch {record["path"]}.')
        if not row.get('files'): errors.append(f'Family row {i}: no source-file receipts.')
    if len(ids)!=len(set(ids)): errors.append('Duplicate family IDs: repeated structures cannot inflate independent n.')
    return {'status':'not_launchable' if errors else 'declared_fields_and_hashes_pass_only',
            'families':len(rows),'errors':errors,
            'scope':'Human source and scientific review is still required; this utility cannot verify eligibility assertions.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('inventory',type=Path);args=ap.parse_args()
    result=audit(json.loads(args.inventory.read_text()),args.inventory.parent)
    print(json.dumps(result,indent=2));raise SystemExit(bool(result['errors']))
