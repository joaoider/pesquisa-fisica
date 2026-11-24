"""T3: I₃ de Rényi-2 invariante, gauge (m=0,J=1) × controle local (V=1, mesmo gap), N=16, protocolo realista."""
import json, time
import numpy as np
from pesquisa_fisica.emergence.ed_quench import evolve, gap, local_interaction_hamiltonian, gap_of, evolve_H
from pesquisa_fisica.emergence.protocolo_t3 import (diag_energies, trotter_trajectory, rho_regions, I3_renyi2_inv_from_rhoABC,
                                                    shadow_batches, batch_purity_fn)

N, m_g, J, V = 16, 0.0, 1.0, 1.0
regs = {"A": [5, 6], "B": [7, 8], "C": [9, 10]}
out = {}
g = gap(N, m_g, J)
lo, hi = 0.0, 6.0
for _ in range(30):
    mid = 0.5 * (lo + hi)
    lo, hi = (mid, hi) if gap_of(local_interaction_hamiltonian(N, 1.0, mid, V)[0]) < g else (lo, mid)
m_c = 0.5 * (lo + hi)
out.update(N=N, gap=g, m_controle=m_c, gap_controle=gap_of(local_interaction_hamiltonian(N, 1.0, m_c, V)[0]))
print(f"gap gauge {g:.4f}; controle V={V} m={m_c:.4f} gap {out['gap_controle']:.4f}", flush=True)

# P1: evolução exata, média t∈[3,6]
times = np.arange(0, 6.0001, 0.5)
stg = evolve(N, m_g, J, times)
Hc, bc = local_interaction_hamiltonian(N, 1.0, m_c, V)
stc = evolve_H(Hc, bc, N, times)
I3g = [I3_renyi2_inv_from_rhoABC(rho_regions(s, N, regs), regs) for s in stg]
I3c = [I3_renyi2_inv_from_rhoABC(rho_regions(s, N, regs), regs) for s in stc]
w = times >= 3
out["P1"] = {"times": times.tolist(), "I3_gauge": I3g, "I3_controle": I3c,
             "avg_gauge": float(np.mean(np.array(I3g)[w])), "avg_controle": float(np.mean(np.array(I3c)[w]))}
print("P1 média t∈[3,6]: gauge %.4f controle %.4f razão %.2f" % (out["P1"]["avg_gauge"], out["P1"]["avg_controle"],
      abs(out["P1"]["avg_gauge"]) / max(abs(out["P1"]["avg_controle"]), 1e-9)), flush=True)
json.dump(out, open("results/t3_protocolo_i3.json", "w"))

# P2: Trotter dt=0,25 até t=4, p=0 e p=0,01 (trajetórias)
dt, steps, R = 0.25, 16, 300
Eg, Ec = diag_energies(N, m_g, J=J), diag_energies(N, m_c, V=V)
rng = np.random.default_rng(2026)
rhos = {}
for p in (0.0, 0.01):
    for name, E in (("gauge", Eg), ("controle", Ec)):
        nr = 1 if p == 0 else R
        rho = sum(np.outer(*(lambda r: (r, r.conj()))(rho_regions(trotter_trajectory(N, E, dt, steps, p, rng), N, regs)[0:0]))
                  for _ in range(0)) if False else None
        acc = np.zeros((64, 64), complex)
        for _ in range(nr):
            psi = trotter_trajectory(N, E, dt, steps, p, rng)
            acc += rho_regions(psi, N, regs)
        rhos[(p, name)] = acc / nr
    Ig, Ic = (I3_renyi2_inv_from_rhoABC(rhos[(p, k)], regs) for k in ("gauge", "controle"))
    out.setdefault("P2", {})[str(p)] = {"I3_gauge": Ig, "I3_controle": Ic, "Delta": Ig - Ic}
    print(f"P2 p={p}: I3 gauge {Ig:.4f} controle {Ic:.4f} Δ {Ig-Ic:.4f}", flush=True)
    json.dump(out, open("results/t3_protocolo_i3.json", "w"))
exact_t4 = {"gauge": I3g[int(np.argmin(abs(times - 4)))], "controle": I3c[int(np.argmin(abs(times - 4)))]}
out["P2"]["exato_t4"] = exact_t4
out["P2"]["retencao"] = out["P2"]["0.01"]["Delta"] / out["P2"]["0.0"]["Delta"]
print("P2 retenção Δ(p=0,01)/Δ(p=0) = %.3f ; exato t=4: Δ=%.4f" % (out["P2"]["retencao"], exact_t4["gauge"] - exact_t4["controle"]), flush=True)
json.dump(out, open("results/t3_protocolo_i3.json", "w"))

# P3: sombras M=10^4, 20 repetições por modelo, p=0,01
M, nb, reps = 10000, 10, 20
est = {}
for name in ("gauge", "controle"):
    t0 = time.time()
    est[name] = [I3_renyi2_inv_from_rhoABC(rhos[(0.01, name)], regs,
                 batch_purity_fn(shadow_batches(rhos[(0.01, name)], 6, M, nb, rng), 6)) for _ in range(reps)]
    print(f"P3 {name}: média {np.mean(est[name]):.4f} desvio {np.std(est[name]):.4f} ({time.time()-t0:.0f}s)", flush=True)
d = np.array(est["gauge"]) - np.array(est["controle"])
sig = float(np.sqrt(np.var(est["gauge"]) + np.var(est["controle"])))
out["P3"] = {"M": M, "reps": reps, "est_gauge": est["gauge"], "est_controle": est["controle"],
             "Delta_medio": float(d.mean()), "sigma": sig, "z": float(abs(d.mean()) / sig)}
print("P3: Δ médio %.4f, σ %.4f, |Δ|/σ = %.2f" % (out["P3"]["Delta_medio"], sig, out["P3"]["z"]), flush=True)
json.dump(out, open("results/t3_protocolo_i3.json", "w"))
