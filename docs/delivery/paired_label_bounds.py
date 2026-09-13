#!/usr/bin/env python3
"""Exact paired P@5 identification bounds; no fitted values or biological labels."""
import itertools

def bounds(pulsar, comparator, labels):
    for prediction in (pulsar, comparator):
        if prediction is not None and (len(prediction) != 5 or len(set(prediction)) != 5):
            raise ValueError('A completed shortlist must contain exactly five unique residues.')
    allowed = {'positive', 'non_contact', 'unknown'}
    if any(value not in allowed for value in labels.values()):
        raise ValueError('Labels must explicitly identify positive, non_contact or unknown.')
    def single(prediction):
        if prediction is None:
            return 0.0, 1.0
        positives = sum(labels.get(r, 'unknown') == 'positive' for r in prediction)
        unknowns = sum(labels.get(r, 'unknown') == 'unknown' for r in prediction)
        return positives / 5, (positives + unknowns) / 5
    if pulsar is None or comparator is None:
        p_lo, p_hi = single(pulsar); c_lo, c_hi = single(comparator)
        return p_lo - c_hi, p_hi - c_lo
    p, c = set(pulsar), set(comparator)
    lo = hi = 0.0
    for r in p | c:
        a = (int(r in p) - int(r in c)) / 5
        label = labels.get(r, 'unknown')
        if label == 'positive':
            lo += a; hi += a
        elif label == 'unknown':
            lo += min(a, 0); hi += max(a, 0)
    return lo, hi

def self_test():
    p=[1,2,3,4,5]; c=[3,4,5,6,7]
    labels={1:'positive',2:'unknown',3:'unknown',4:'non_contact',5:'positive',6:'unknown',7:'non_contact'}
    lo,hi=bounds(p,c,labels)
    possible=[]
    unknown=[r for r,v in labels.items() if v=='unknown']
    for assigned in itertools.product(['positive','non_contact'],repeat=len(unknown)):
        completed=labels | dict(zip(unknown,assigned))
        value=sum(completed[r]=='positive' for r in p)/5-sum(completed[r]=='positive' for r in c)/5
        possible.append(value)
    assert abs(lo-min(possible))<1e-12 and abs(hi-max(possible))<1e-12
    assert bounds(p,p,{})==(0.0,0.0), 'Shared unknowns must cancel.'
    assert bounds(p,c,{})==(-.4,.4), 'Absent labels remain unknown.'
    assert bounds(None,c,{})==(-1.0,1.0), 'Failed method cannot become a zero-accuracy comparator.'
    assert bounds(p,None,{r:'positive' for r in p})==(0.0,1.0)
    for bad in ([1,2,3,4,4], [1,2,3,4]):
        try: bounds(bad,c,{})
        except ValueError: pass
        else: raise AssertionError('Invalid shortlist was accepted.')
    print('PASS: exhaustive shared-label bounds; missing labels stay unknown; failed-method bounds; fixed-five validity.')

if __name__=='__main__': self_test()
