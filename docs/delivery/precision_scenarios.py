#!/usr/bin/env python3
"""Analytical study-design scenarios only; no protein observations are generated."""
import csv, hashlib, json, math, platform
from pathlib import Path
import scipy
from scipy.integrate import quad
from scipy.stats import chi2, nct, norm, t

ROOT = Path(__file__).resolve().parent
NS = (12, 20, 30, 40, 60, 100)
SDS = (0.15, 0.25, 0.40, 0.60)
DELTAS = (0.0, 0.05, 0.10, 0.20)
ALPHA = 0.05
POINT_TARGET = 0.10

def half_width(n, sd):
    return float(t.ppf(1 - ALPHA / 2, n - 1) * sd / math.sqrt(n))

def joint_success_probability(n, sd, delta):
    # Under independent normal family differences, sample mean and sample
    # variance are independent. Integrate the mean threshold over chi-square
    # sample-variance uncertainty. This is NOT a biological power estimate.
    df = n - 1
    critical = float(t.ppf(1 - ALPHA / 2, df))
    se = sd / math.sqrt(n)
    def integrand(v):
        threshold = max(POINT_TARGET, critical * se * math.sqrt(v / df))
        return norm.sf((threshold - delta) / se) * chi2.pdf(v, df)
    value, numerical_error = quad(integrand, 0.0, math.inf,
                                  epsabs=1e-10, epsrel=1e-10, limit=250)
    return float(value), float(numerical_error)

def main():
    rows = []
    for n in NS:
        for sd in SDS:
            for delta in DELTAS:
                joint, integration_error = joint_success_probability(n, sd, delta)
                positive = float(nct.sf(t.ppf(.975, n - 1), n - 1,
                                        delta * math.sqrt(n) / sd))
                rows.append(dict(independent_families=n, assumed_paired_sd=sd,
                    assumed_mean_gain=delta, nominal_t_half_width=half_width(n, sd),
                    probability_positive_95percent_interval=positive,
                    probability_positive_interval_and_estimate_at_least_0p10=joint,
                    integration_absolute_error_estimate=integration_error))
    with (ROOT / 'precision-scenarios.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, rows[0].keys()); writer.writeheader(); writer.writerows(rows)
    requirements=[]
    for sd in SDS:
        for width in (.10, .05):
            n=next(n for n in range(3, 10001) if half_width(n, sd) <= width)
            requirements.append(dict(assumed_paired_sd=sd, desired_nominal_half_width=width,
                                     required_independent_families=n))
    checks = {
        '96_scenarios': len(rows) == 96,
        'zero_gain_positive_interval_probability_0p025': all(abs(r['probability_positive_95percent_interval']-.025)<1e-10 for r in rows if r['assumed_mean_gain']==0),
        'true_target_joint_success_no_more_than_one_half': all(r['probability_positive_interval_and_estimate_at_least_0p10']<=.5+1e-10 for r in rows if r['assumed_mean_gain']==.10),
        'joint_probability_not_greater_than_positive_interval_probability': all(r['probability_positive_interval_and_estimate_at_least_0p10']<=r['probability_positive_95percent_interval']+1e-9 for r in rows),
        'integration_error_below_1e_minus_8': max(r['integration_absolute_error_estimate'] for r in rows)<1e-8,
        'width_decreases_with_n': all(half_width(20, s)>half_width(40,s)>half_width(100,s) for s in SDS),
    }
    assert all(checks.values()), checks
    result = dict(status='analytical_planning_scenarios_not_biological_evidence',
        assumptions=['Independent, identically distributed approximately normal family-level paired differences for t/noncentral-t planning.',
                     'All sigma and delta values are assumed; no independent empirical sigma is available.',
                     'Precision planning is not a guarantee of realized interval width or coverage for discrete/skewed family outcomes.'],
        formula='half_width = t_(0.975,n-1) * assumed_paired_sd / sqrt(n)',
        required_n_for_nominal_precision=requirements,
        bounded_difference_sensitivity=dict(range=[-1,1],
            formula='sqrt(2*log(2/alpha)/n)',
            n_for_half_width_0p10=math.ceil(2*math.log(2/ALPHA)/.1**2),
            caveat='Hoeffding bound assumes independent bounded family differences; it does not establish representative sampling.'),
        tests=checks,python=platform.python_version(),scipy=scipy.__version__,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        csv_sha256=hashlib.sha256((ROOT/'precision-scenarios.csv').read_bytes()).hexdigest())
    (ROOT/'precision-summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'tests':checks,'required_n':requirements,'scenario_count':len(rows)},indent=2))

if __name__ == '__main__': main()
