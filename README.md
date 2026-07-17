# Entanglement-Geometric Signatures of Hyper-Stiff Cosmological Phases

[![arXiv](https://shields.io)](https://arxiv.org)
[![License: MIT](https://shields.io)](https://opensource.org)

Official repository for the paper **"Entanglement-Geometric Signatures of Hyper-Stiff Cosmological Phases: A Unified Framework Linking Brane Rupture, Confinement Energy, and Mutual-Information Geometry"** by Paul Jarvis.

This repository provides the numerical suite used to demonstrate how early-universe cosmological equations of state map directly onto distinct geometric regimes reconstructed from pure quantum entanglement correlations.

---

## Framework Overview

This framework establishes a model-independent dictionary between early-universe dynamics and quantum information metrics by mapping cosmological fluid behavior onto corresponding spin-chain analogue ground states:

*   **Radiation-like Phase ($w = 1/3$):** Modeled via a critical $XX$ chain ($J_z = 0$), yielding a smooth, extended, and gapless emergent geometry.
*   **Matter-like Phase ($w = 0$):** Modeled via a weakly gapped $XXZ$ chain ($J_z = 0.5$), retaining spatial extension with an emergent mass scale.
*   **Hyper-Stiff Rupture Phase ($w \gg 1$):** Modeled via a strongly gapped Ising-like chain ($J_z = 2.0$), corresponding to a highly connected, ultra-dense, and collapsing geometric regime.

---

## Key Physical Results

1. **The Inverse Scaling Law Resolved:** While radiation and matter analogues show standard geometric scaling ($\lambda_1 \propto 1/N^2$), hyper-stiff phases show a monotonic *increase* in $\lambda_1$ with system size. Because the mutual information decay flatlines ($\alpha \approx 0$), the underlying state approaches a maximally connected graph whose Laplacian naturally yields $\lambda_1 \propto N$, signaling a physical transition to an ultra-dense, collapsing information topology.
2. **Clear Phase Clustering:** Plotting the spectral gap against the decay exponent ($\lambda_1$ vs. $\alpha$) or localization (IPR vs. $\alpha$) isolates three non-overlapping phase clusters, consistent with entanglement structure tracking cosmic epoch.
3. **Noise and Stability Controls:** Bootstrap resampling and random-matrix null models confirm that the spectral gap signature of the hyper-stiff phase is not replicated by statistical noise.

---

## Repository Structure

* `/code/testing_code.py`: Production script containing the exact diagonalization routines, Laplacian solvers, bootstrap resamplers, and null-model controls. Generates Figures 1-10 into `/figures`.
* `/figures`: Generated publication-grade plots (Figures 1-10; Figure 9 is not referenced in the manuscript and is superseded by Figure 10's corrected version of the same diagnostic).
* `/main.tex`: LaTeX manuscript source.
* `/main.pdf`: The compiled research paper.

---

## Quick Start & Reproducibility

To re-run the pipeline and reproduce the figures from the manuscript:

```bash
# Clone the repository
git clone https://github.com/ZaPpi3/hyperstiff-entanglement-geometry
cd hyperstiff-entanglement-geometry

# Install requirements
pip install numpy scipy matplotlib qutip

# Run the complete pipeline (writes into figures/)
python code/testing_code.py
```

Requires `qutip` for exact diagonalization of the spin-chain Hamiltonians; the largest system size used ($N=14$) is not fast.

To rebuild the manuscript PDF from source:
```bash
tectonic main.tex
```
