import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from scipy.linalg import eigh
from qutip import *

# ============================================================
# CONFIGURATION
# ============================================================

Ns = (8, 10, 12, 14)
PHASES = ("radiation_like", "matter_like", "hyperstiff_like")
BOOTSTRAP_SAMPLES = 200

# ============================================================
# CORE UTILITIES
# ============================================================

def op_on_site(op, i, N):
    ops = [qeye(2)] * N
    ops[i] = op
    return tensor(ops)

def mutual_information_matrix(psi0, N):
    cache = {}
    def rho1(i):
        if i not in cache:
            cache[i] = psi0.ptrace(i)
        return cache[i]

    I = np.zeros((N, N))
    for i, j in product(range(N), range(N)):
        if i == j:
            continue
        if I[j, i] != 0:
            I[i, j] = I[j, i]
            continue
        rho_i = rho1(i)
        rho_j = rho1(j)
        rho_ij = psi0.ptrace([i, j])
        I[i, j] = entropy_vn(rho_i) + entropy_vn(rho_j) - entropy_vn(rho_ij)
    return I

def laplacian_from_MI(I):
    A = I / np.max(I)
    np.fill_diagonal(A, 0)
    D = np.diag(np.sum(A, axis=1))
    return D - A

def inverse_participation_ratio(v):
    v2 = v**2
    return np.sum(v2**2) / (np.sum(v2)**2)

# ============================================================
# TOY HAMILTONIANS
# ============================================================

def build_hamiltonian(N, phase):
    sx, sy, sz = sigmax(), sigmay(), sigmaz()

    if phase == "radiation_like":
        Jx, Jy, Jz = 1.0, 1.0, 0.0
    elif phase == "matter_like":
        Jx, Jy, Jz = 1.0, 1.0, 0.5
    elif phase == "hyperstiff_like":
        Jx, Jy, Jz = 0.2, 0.2, 2.0

    H = 0
    for i in range(N - 1):
        H += (
            Jx * op_on_site(sx, i, N) * op_on_site(sx, i+1, N) +
            Jy * op_on_site(sy, i, N) * op_on_site(sy, i+1, N) +
            Jz * op_on_site(sz, i, N) * op_on_site(sz, i+1, N)
        )
    return H

def ground_state(N, phase):
    H = build_hamiltonian(N, phase)
    evals, evecs = H.eigenstates(eigvals=1, sparse=True)
    return evecs[0], evals[0]

# ============================================================
# ANALYSIS
# ============================================================

def MI_decay_exponent(I, N):
    rs = np.arange(1, N)
    Ivals, rlist = [], []
    for r in rs:
        vals = [I[i, i+r] for i in range(N-r)]
        vals = [v for v in vals if v > 0]
        if len(vals) == 0:
            continue
        Ivals.extend(vals)
        rlist.extend([r]*len(vals))
    log_r = np.log(rlist)
    log_I = np.log(Ivals)
    alpha, _ = np.polyfit(log_r, log_I, 1)
    return alpha

def analyze_phase(N, phase):
    psi0, _ = ground_state(N, phase)
    I = mutual_information_matrix(psi0, N)
    L = laplacian_from_MI(I)
    evals, vecs = eigh(L)
    lam1 = evals[1]
    alpha = MI_decay_exponent(I, N)
    ipr = inverse_participation_ratio(vecs[:, 1])
    return lam1, alpha, ipr, I, evals, vecs

def bootstrap_lambda1(I, samples=200):
    lam1s = []
    for _ in range(samples):
        noise = np.random.normal(0, 0.05 * np.mean(I[I > 0]), size=I.shape)
        I_noisy = np.abs(I + noise)
        L = laplacian_from_MI(I_noisy)
        evals, _ = eigh(L)
        lam1s.append(evals[1])
    return np.mean(lam1s), np.std(lam1s)

def random_null_model(N):
    R = np.random.rand(N, N)
    R = (R + R.T) / 2
    np.fill_diagonal(R, 0)
    L = laplacian_from_MI(R)
    evals, _ = eigh(L)
    return evals[1]

# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    results = []

    for N in Ns:
        for phase in PHASES:
            lam1, alpha, ipr, I, evals, vecs = analyze_phase(N, phase)
            lam1_bs_mean, lam1_bs_std = bootstrap_lambda1(I, BOOTSTRAP_SAMPLES)
            lam1_null = random_null_model(N)

            results.append({
                "N": N,
                "phase": phase,
                "lam1": lam1,
                "lam1_bs_mean": lam1_bs_mean,
                "lam1_bs_std": lam1_bs_std,
                "lam1_null": lam1_null,
                "alpha": alpha,
                "ipr": ipr,
            })

    # ============================================================
    # 10 PUBLICATION-READY FIGURES
    # ============================================================

    # 1. MI decay curves
    plt.figure()
    for phase in PHASES:
        for N in Ns:
            r = next(r for r in results if r["N"] == N and r["phase"] == phase)
            I = mutual_information_matrix(ground_state(N, phase)[0], N)
            rs = np.arange(1, N)
            mean_I = [np.mean([I[i, i+d] for i in range(N-d)]) for d in rs]
            plt.loglog(rs, mean_I, label=f"{phase}, N={N}")
    plt.xlabel("r")
    plt.ylabel("⟨I(i,i+r)⟩")
    plt.title("MI decay curves")
    plt.legend()
    plt.savefig("FIG1_MI_decay.png", dpi=300)
    plt.close()

    # 2. Eigenmode comparison
    plt.figure()
    for phase in PHASES:
        r = next(r for r in results if r["N"] == 14 and r["phase"] == phase)
        _, _, _, _, evals, vecs = analyze_phase(14, phase)
        plt.plot(vecs[:, 1], label=phase)
    plt.title("First nonzero eigenmode comparison (N=14)")
    plt.legend()
    plt.savefig("FIG2_eigenmodes.png", dpi=300)
    plt.close()

    # 3. Spectral density comparison
    plt.figure()
    for phase in PHASES:
        _, _, _, _, evals, _ = analyze_phase(14, phase)
        plt.plot(evals, 'o-', label=phase)
    plt.title("Spectral density comparison (N=14)")
    plt.legend()
    plt.savefig("FIG3_spectrum.png", dpi=300)
    plt.close()

    # 4. Finite-size scaling
    plt.figure()
    for phase in PHASES:
        xs = [1/r["N"]**2 for r in results if r["phase"] == phase]
        ys = [r["lam1"] for r in results if r["phase"] == phase]
        plt.plot(xs, ys, 'o-', label=phase)
    plt.xlabel("1/N²")
    plt.ylabel("λ₁")
    plt.title("Finite-size scaling of λ₁")
    plt.legend()
    plt.savefig("FIG4_finite_size.png", dpi=300)
    plt.close()

    # 5. Phase diagram λ1 vs α
    plt.figure()
    for phase in PHASES:
        xs = [r["alpha"] for r in results if r["phase"] == phase]
        ys = [r["lam1"] for r in results if r["phase"] == phase]
        plt.scatter(xs, ys, label=phase)
    plt.xlabel("MI decay exponent α")
    plt.ylabel("λ₁")
    plt.title("Phase diagram: λ₁ vs α")
    plt.legend()
    plt.savefig("FIG5_phase_diagram.png", dpi=300)
    plt.close()

    # 6. Bootstrap stability
    plt.figure()
    for phase in PHASES:
        Ns_plot = [r["N"] for r in results if r["phase"] == phase]
        lam1_plot = [r["lam1_bs_mean"] for r in results if r["phase"] == phase]
        err_plot = [r["lam1_bs_std"] for r in results if r["phase"] == phase]
        plt.errorbar(Ns_plot, lam1_plot, yerr=err_plot, fmt='o-', label=phase)
    plt.xlabel("N")
    plt.ylabel("λ₁ (bootstrap mean ± σ)")
    plt.title("Bootstrap stability of λ₁")
    plt.legend()
    plt.savefig("FIG6_bootstrap.png", dpi=300)
    plt.close()

    # 7. Null model comparison
    plt.figure()
    for phase in PHASES:
        Ns_plot = [r["N"] for r in results if r["phase"] == phase]
        lam1_plot = [r["lam1"] for r in results if r["phase"] == phase]
        null_plot = [r["lam1_null"] for r in results if r["phase"] == phase]
        plt.plot(Ns_plot, lam1_plot, 'o-', label=f"{phase} λ₁")
        plt.plot(Ns_plot, null_plot, '--', label=f"{phase} null")
    plt.xlabel("N")
    plt.ylabel("λ₁")
    plt.title("Null-model comparison")
    plt.legend()
    plt.savefig("FIG7_null_model.png", dpi=300)
    plt.close()

    # 8. IPR vs N
    plt.figure()
    for phase in PHASES:
        Ns_plot = [r["N"] for r in results if r["phase"] == phase]
        ipr_plot = [r["ipr"] for r in results if r["phase"] == phase]
        plt.plot(Ns_plot, ipr_plot, 'o-', label=phase)
    plt.xlabel("N")
    plt.ylabel("IPR")
    plt.title("Localization (IPR) vs N")
    plt.legend()
    plt.savefig("FIG8_IPR_vs_N.png", dpi=300)
    plt.close()

    # 9. IPR vs λ1
    plt.figure()
    for phase in PHASES:
        xs = [r["lam1"] for r in results if r["phase"] == phase]
        ys = [r["ipr"] for r in results if r["phase"] == phase]
        plt.scatter(xs, ys, label=phase)
    plt.xlabel("λ₁")
    plt.ylabel("IPR")
    plt.title("Localization vs mass scale")
    plt.legend()
    plt.savefig("FIG9_IPR_vs_lam1.png", dpi=300)
    plt.close()


    # 10. IPR vs alpha (corrected diagnostic)
    plt.figure()
    for phase in PHASES:
        xs = [r["alpha"] for r in results if r["phase"] == phase]
        ys = [r["ipr"] for r in results if r["phase"] == phase]
        plt.scatter(xs, ys, label=phase)
    plt.xlabel("MI decay exponent alpha")
    plt.ylabel("IPR")
    plt.title("Localization vs MI decay exponent alpha")
    plt.legend()
    plt.savefig("FIG10_alpha_vs_lam1.png", dpi=300)
    plt.close()
