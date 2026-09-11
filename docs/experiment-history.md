# Experiment history

The initial fixed synthetic model used κ=10, β=1 and two coordinates. Its strain diagnostics prompted an exploratory stiffness sweep (10, 30, 100, 300, 1000), domain checks and coordinate refinement. κ=100 was selected for the main illustration after inspecting strain, not biological accuracy. Every recorded setting is retained. These are exploratory verification experiments, not preregistered biological tests.

Initial macOS NumPy matrix multiplications emitted floating-point warnings even when outputs were finite. Independent scalar-energy and spectral calculations reproduced the results. The released source explicitly uses inspectable einsum contractions and an independent small-matrix exponential reference. It does not suppress warnings or rewrite source code at runtime. The original mathematical equations are unchanged; released code and all scientific outputs are retested in a clean environment.

Floating-point values are compared with declared tolerances. Runtime, timestamps and platform metadata are not expected to match bit for bit. There are no random simulation trajectories or stochastic biological replicates in these tests.
