\# Entanglement-Geometric Signatures of Hyper-Stiff Cosmological Phases



\[!\[arXiv](https://shields.io)](https://arxiv.org)

\[!\[License: MIT](https://shields.io)](https://opensource.org)



Official repository for the paper \*\*“Entanglement-Geometric Signatures of Hyper-Stiff Cosmological Phases: A Unified Framework Linking Brane Rupture, Confinement Energy, and Mutual-Information Geometry” (2026)\*\* by Paul Jarvis \[1.1].



This repository provides the numerical suite used to demonstrate how early-universe cosmological equations of state map directly onto distinct geometric regimes reconstructed from pure quantum entanglement correlations \[1.1].



\---



\## 🔍 Framework Overview



This framework establishes a precise, model-independent dictionary between early-universe dynamics and quantum information metrics by mapping cosmological fluid behavior onto corresponding spin-chain analogue ground states \[1.1, 1.2]:



\*   \*\*Radiation-like Phase ($w = 1/3$):\*\* Modeled via a critical $XX$ chain ($J\_z = 0$), yielding a smooth, extended, and gapless emergent geometry \[1.1, 1.2, 1.4].

\*   \*\*Matter-like Phase ($w = 0$):\*\* Modeled via a weakly gapped $XXZ$ chain ($J\_z = 0.5$), retaining spatial extension with an emergent mass scale \[1.1, 1.2, 1.4].

\*   \*\*Hyper-Stiff Rupture Phase ($w \\gg 1$):\*\* Modeled via a strongly gapped Ising-like chain ($J\_z = 2.0$), corresponding to a highly connected, ultra-dense, and collapsing geometric regime \[1.1, 1.2, 1.4].



\---



\## 🧠 Key Physical Results



1\. \*\*The Inverse Scaling Law Resolved:\*\* While radiation and matter analogues show standard geometric scaling ($\\lambda\_1 \\propto 1/N^2$), hyper-stiff phases show a monotonic \*increase\* in $\\lambda\_1$ with system size \[1.10, 1.11]. Because the mutual information decay flatlines ($\\alpha \\approx 0$), the underlying state approaches a maximally connected graph whose Laplacian naturally yields $\\lambda\_1 \\propto N$, signaling a physical transition to an ultra-dense, collapsing information topology \[1.11].

2\. \*\*Clear Phase Clustering:\*\* Plotting the spectral gap against the decay exponent ($\\lambda\_1$ vs. $\\alpha$) or localization ($\\text{IPR}$ vs. $\\alpha$) isolates three unassailable, non-overlapping phase clusters, confirming that entanglement structures uniquely classify cosmic epochs \[1.11, 1.12].

3\. \*\*Noise and Stability Controls:\*\* Bootstrap resampling and random-matrix null models confirm that the massive spectral gap signature of the hyper-stiff phase cannot be replicated by statistical noise \[1.11].



\---



\## 📁 Repository Structure



\*   `/code` : Production script (`main.py`) containing the exact diagonalization routines, Laplacian solvers, bootstrap resamplers, and null-model controls \[1.10].

\*   `/output` : Generated publication-grade plots (Figures 1–9) \[1.10].

\*   `manuscript.pdf` : The compiled research paper \[1.1].



\---



\## 🚀 Quick Start \& Reproducibility



To re-run the pipeline and completely reproduce Figures 1 through 9 from the manuscript \[1.10]:



```bash

\# Clone the repository

git clone https://github.com

cd Emergent-Geometry-from-Entanglement



\# Install requirements

pip install numpy scipy matplotlib



\# Run the complete pipeline

python code/main.py

```



\---



\## 📜 Citation



```bibtex

@article{jarvis2026cosmoentanglement,

&#x20; title={Entanglement-Geometric Signatures of Hyper-Stiff Cosmological Phases: A Unified Framework Linking Brane Rupture, Confinement Energy, and Mutual-Information Geometry},

&#x20; author={Jarvis, Paul},

&#x20; journal={arXiv preprint arXiv:2606.XXXXX},

&#x20; year={2026},

&#x20; month={June}

}

```



