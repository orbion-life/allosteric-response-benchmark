# Protein response validation and quantum execution evidence

14 September 2026. This review supports the revised Pulsar proposal. It distinguishes completed measurements, literature support, calculated execution budgets and proposed work. It does not establish selection by the jury, a biological advantage or practical quantum advantage.

## Decision supported by the new work

The strongest defensible proposal is a fourteen-week investigation with a working four-target classical delivery baseline and a newly constructed protein-derived quantum estimator. The exact spectral representation closes the previous missing-circuit-construction evidence for KRAS. It does not make complete protein quantum execution practical: the same spectral preprocessing already makes the classical calculation fast, and the evaluated quantum schedule remains large.

Biological usefulness remains an empirical question. The completed KRAS and GRB2 comparisons do not establish an improvement over structural controls. An additional distinct-family extension has a frozen protocol and input-qualification records but no completed prediction or outcome comparison. It supplies no new biological result and is not counted as a failed model test.

The [official Phase 1 assessment](https://quantumaiportal.thequantuminsider.com/wp-content/uploads/2026/04/2026-04-06-Assessment-Criteria-VF.pdf) evaluates relevance and impact, technical approach and innovation, feasibility, validation, hybrid integration and team capability. A credible plan matters at this phase; already established quantum advantage is not an explicit prerequisite. Nevertheless, an inexpensive classical solution limits the rationale for presenting the spectral estimator itself as a route to speed advantage.

## The claim has several distinct evidential steps

| Claim | Evidence now available | Defensible conclusion |
|---|---|---|
| Static structures can produce complete response matrices at protein size | Four Gaussian calculations, up to 694 residues; largest recorded run 249.19 seconds and 988 MB | Computation is demonstrated within this Gaussian approximation. |
| The reduced representation preserves the declared outputs | Twelve target/time queries pass 0.002; maximum error 0.000902 | Preservation is supported at the three declared times. |
| It preserves general transient responses | Expanded step and pulse comparisons fail | Arbitrary-time response preservation remains unsupported. |
| The model preserves original nonlinear physics | Gaussian triangle discrepancy 0.4093 | Original-model fidelity is not established by the Gaussian result. |
| The response can be represented by a protein-derived quantum circuit | All KRAS observable preparations and filters constructed; representative 20-qubit Grover simulation agrees with a fresh matrix-exponential reference | Circuit construction and selected ideal semantics are demonstrated. |
| It improves biological prioritization | Gaussian correlations are positive in two completed assays, but structural controls remain competitive | Incremental biological benefit is not established. |
| Complete quantum execution is practical or faster | Constructive costs remain large; classical spectral calculation is fast | Neither practical full quantum execution nor speed advantage is established. |

Prior evidence is retained in [v0.7.0](https://github.com/orbion-life/allosteric-response-benchmark/releases/tag/v0.7.0) and [v0.6.0](https://github.com/orbion-life/allosteric-response-benchmark/releases/tag/v0.6.0). Changing the estimator does not repair a physical or biological approximation upstream.

## What the literature supports

[Faure et al. (2022)](https://www.nature.com/articles/s41586-022-04586-4) provide binding and abundance measurements interpreted through thermodynamic models. Their work supports separating binding effects from folding effects when constructing a functional endpoint. It does not imply that a static contact-response score predicts those effects.

[Weng et al. (2024)](https://www.nature.com/articles/s41586-023-06954-0) supply a relevant KRAS energetic landscape. The existing Pulsar comparison uses a partner-specific binding endpoint and retains the structure–assay nucleotide-state mismatch. An inferred binding-energy magnitude, a known ligand contact and direct regulatory activity are different outcomes; they must not be pooled into one accuracy claim.

[Martí-Aranda and Lehner (2026)](https://www.nature.com/articles/s41467-026-71005-x) report seven interaction maps across five homologous PDZ domains, with both distance dependence and protein-specific distal effects. This strengthens the case for a distance control and for testing transfer. Seven interactions are not seven independent protein families, and one map reuses the earlier PSD95-PDZ3 data. The paper's displayed structures include AlphaFold predictions, so its existence does not establish admission under Pulsar's current exact experimental apo-structure protocol.

[Mohtashim, Sajjan and Kais (2026)](https://pubs.acs.org/doi/abs/10.1021/jacs.6c08053) establish relevant quantum-walk protein-network prior art, including a reported hardware proof of principle. Their centrality observable differs from Pulsar's signed contact-energy response. Their abstract does not establish a hardware cost or predictive advantage for Pulsar, and the 150-protein structural comparison must not be presented as 150 independent allosteric-function validations.

[Grinko et al. (2021)](https://arxiv.org/abs/1912.05559) motivate amplitude estimation without quantum phase estimation. Their query-complexity result concerns an oracle model and comparison with Monte Carlo sampling. Pulsar separately implements a fresh-batch confidence procedure and charges actual preparations, inverses and reflections. A lower sampling query count does not establish a speedup over deterministic spectral matrix multiplication.

[Zhang and Yuan (2024)](https://www.nature.com/articles/s41534-024-00835-8) make data-access circuit costs relevant to the comparison. [Litinski (2019)](https://quantum-journal.org/papers/q-2019-03-05-128/) shows why logical circuits must be mapped to an error-corrected architecture with space–time and magic-state costs before a physical resource claim is complete. Neither paper supplies a free hardware implementation for this proposal. [Hoefler, Häner and Troyer (2023)](https://arxiv.org/abs/2307.00523) provide a broader practical-advantage analysis; it is context rather than a theorem ruling out every future Pulsar representation.

## The new quantum calculation

The spectral construction diagonalizes the saved positive operator and transforms every observable. Controlled signed state preparation and an attenuation flag encode the complete reduced kernel. The response is recovered by rescaling a Hadamard probability and subtracting the static covariance. Every eigenmode is retained. Both the classical and quantum branches pay for diagonalization and observable transformation.

| Quantity | Rank-nine triangle | Rank-347 KRAS |
|---|---:|---:|
| Active circuit qubits | 6 | 11 |
| Qubits including Grover reflection helpers | 10 | 20 |
| Native CX gates per complete preparation | 92 | 3,068 |
| All pair/time outputs | 18 | 41,583 |
| Realized ideal binomial readouts | 20,992 | 63,497,472 |
| Constructive campaign T count | 141,682,054,144 | 56,450,203,649,724,416 |
| Constructive campaign CX count | 118,141,952 | 32,712,522,454,784 |
| T ratio to fixed sufficient spectral sampling plan | 0.0372844 | 0.0133605 |
| Maximum Grover power | 212 | 633 |
| Maximum CX count in one coherent shot | 44,400 | 3,921,971 |
| Complete response error allowance | 0.00145969 | 0.00196847 |

These counts include all unordered pairs, including diagonals, at three times; symmetry reconstructs the full matrices. Every ideal interval covered its reference in the recorded realization. The 95% joint confidence argument uses conditional exact binomial intervals and failure allocations that sum to at most 0.05 across queries and rounds. Empirical coverage of one realization is a diagnostic, not the confidence proof.

The adaptive likelihood simulations are classical binomial draws from exact ideal probabilities. The scheduling routine receives sampled counts rather than the known amplitude. Complete T counts compose actual synthesized binary rotation banks with stored schedules. This is a constructive, unoptimized implementation, not a lower bound, a globally optimal quantum algorithm or an executed fault-tolerant campaign. The fixed Hoeffding comparator is a sufficient sampling budget, not the cheapest conceivable sampling procedure. Approximate compiled gates can change the realized adaptive schedule slightly; the stated totals are not guaranteed physical-runtime bounds.

On the recorded M4 Pro calculation, KRAS eigendecomposition, all observable transformations and three complete response matrices took 13.73 ms once the reduced operator existed. Matrix evaluation alone took 0.808 ms. Direct matrix exponentiation agreed within 1.19 × 10⁻¹³. Model fitting and projection remain separate costs. The spectral preprocessing therefore already solves the deterministic reduced problem efficiently. A future speed-advantage claim needs a different evidence-supported computational bottleneck, rather than comparison only with an older expensive quantum encoding.

A complete physical estimate still needs a selected device or code, routed circuits, gate durations and errors, logical error allocation, distillation factories, classical feedback latency, batching, queue/access limits and the complete output set. The revised proposal commits to qualifying a smaller hardware demonstration first. It does not promise that a 63-million-readout full KRAS campaign will run during the sprint.

## Biological usefulness requires a better comparative result

The completed GRB2 evaluation gives Spearman correlations 0.381 for Gaussian response, 0.335 for harmonic response and 0.481 for contact degree on 21 common distal sites. Gaussian minus harmonic is about 0.045, with a paired 95% interval spanning zero; Gaussian and harmonic produce the same top five. The available evidence therefore supports a reproducible comparison, not a demonstrated prioritization gain. The KRAS comparison similarly retains a stronger degree correlation than the Gaussian score.

The appropriate decision is to keep the model and controls fixed for an eligible external evaluation and to report every admitted result. Structural-model changes, endpoint choices, time selection and newly inspected labels must be recorded as development exposure. A missing pocket annotation is unknown, not a demonstrated nonfunctional pocket. A nonsignificant residue effect also does not establish functional equivalence. Supported negative-pocket evidence remains a separate requirement of the challenge's specificity test.

The [statistical contract audit](validation-contract-statistics.md) corrects an additional weakness. Four independent nonzero paired target differences have a minimum conventional exact one-sided sign probability of 0.0625; two have a floor of 0.25. More mutations or alternative structures do not increase the number of independent targets. These facts were independently verified by enumerating all precision@5 difference vectors for two and four targets. They do not rule out every parametric method, but such methods need their own justified assumptions.

The revised contract therefore separates the fixed four-target deliverable from an externally powered generalization study. A mean precision@5 improvement of 0.10 remains a practical target; it is not promised significance. The external cohort must be sized for variation, dependence, label coverage, multiplicity and a stated power or precision objective before its outcomes are inspected. Five or six targets merely make certain exact p-values attainable; they are not automatically adequate study sizes. Unknown selected labels require identification bounds, with shared candidates cancelling in paired comparisons.

## Selection strategy and submission decisions

Lead with the delivered protein-scale calculations and the constructed KRAS circuit. Explain the input, perturbation, signed output and five-residue decision before introducing gate counts. Present the Gaussian approximation explicitly and place its failed nonlinear and transient checks beside the corresponding claims. Promise the four matrices, four shortlists, maps, matched controls, reproducible software and a decision record within fourteen weeks. Treat biological gain and hardware value as measured acceptance decisions.

The proposed team is **Orbion Pulsar**, led by Aniruddh Goteti at aniruddh.goteti@orbion.life. The recommended weekly allocations are 16 hours for Aniruddh, 10 for Çağlar, 8 for Promit and 4 for Dennis, totaling 532 person-hours over fourteen weeks. These are proposed responsibilities and resource decisions, not confirmed availability or institutional affiliations.

The administrative review distinguishes official profile fields from recommended planning details. Internal draft declarations are supplied separately to the team. Exact affiliations, factual eligibility, disclosure rights and contributor agreement remain to be confirmed. Draft declarations are ready for review. Neither delegated planning authority nor agreement with the proposal establishes those existing facts; challenge terms acceptance and actual submission remain separate actions.

The new evidence makes the implementation plan more concrete. A 10/10 or a claim that both scientific limits are resolved would still exceed the evidence.
