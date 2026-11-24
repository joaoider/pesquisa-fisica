"""v0.22: quench exato N=20 — I(L:R), I₃ e S de região única nas álgebras qubit e inv; gauge vs livre (mesmo gap)."""
import json, time
import numpy as np
from pesquisa_fisica.emergence.ed_quench import evolve, gap
from pesquisa_fisica.emergence.gauge_algebra import algebra_entropy, algebra_I3
from pesquisa_fisica.emergence.schwinger_dmrg import free_gap

N = 20
times = np.round(np.arange(0, 6.0001, 0.5), 3)
L, M, R = list(range(0, 6)), list(range(6, 14)), list(range(14, 20))
A, B, C = [7, 8], [9, 10], [11, 12]
half = list(range(10))


def match(target, lo=0.0, hi=6.0):
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if free_gap(N, mid) < target else (lo, mid)
    return 0.5 * (lo + hi)


def measure(st):
    rec = {k: [] for k in ("S_half_q", "S_half_i", "S_M_q", "S_M_i", "ILR_q", "ILR_i", "I3_q", "I3_i")}
    for s in st:
        rec["S_half_q"].append(algebra_entropy(s, N, [half], (), "qubit"))
        rec["S_half_i"].append(algebra_entropy(s, N, [half], (), "inv"))
        rec["S_M_q"].append(algebra_entropy(s, N, [M], (), "qubit"))
        rec["S_M_i"].append(algebra_entropy(s, N, [M], (), "inv"))
        SL, SR = algebra_entropy(s, N, [L], (), "inv"), algebra_entropy(s, N, [R], (), "inv")
        rec["ILR_q"].append(SL + SR - rec["S_M_q"][-1])                    # pureza: S(L∪R) = S(M) na álgebra completa
        rec["ILR_i"].append(SL + SR - algebra_entropy(s, N, [L, R], (), "inv"))
        rec["I3_q"].append(algebra_I3(s, N, A, B, C, "qubit"))
        rec["I3_i"].append(algebra_I3(s, N, A, B, C, "inv"))
    return rec


out = []
for m, J in [(0.0, 1.0), (0.5, 2.0)]:
    t0 = time.time()
    g = gap(N, m, J); mf = match(g)
    row = {"m": m, "J": J, "gap": g, "m_free": mf, "times": times.tolist(),
           "gauge": measure(evolve(N, m, J, times)), "free": measure(evolve(N, mf, 0.0, times))}
    row["secs"] = round(time.time() - t0, 1)
    out.append(row)
    json.dump(out, open("results/quench_invariante_N20.json", "w"))
    G, F = row["gauge"], row["free"]
    print(f"m={m} J={J} gap={g:.3f} ({row['secs']}s)")
    print("  max|S_q−S_i| meia/bloco:", max(abs(a - b) for a, b in zip(G['S_half_q'] + F['S_half_q'], G['S_half_i'] + F['S_half_i'])),
          max(abs(a - b) for a, b in zip(G['S_M_q'] + F['S_M_q'], G['S_M_i'] + F['S_M_i'])))
    print(f"  t=6  I(L:R) qubit g/f = {G['ILR_q'][-1]:.4f}/{F['ILR_q'][-1]:.4f} = {G['ILR_q'][-1]/F['ILR_q'][-1]:.3f}"
          f" | inv g/f = {G['ILR_i'][-1]:.4f}/{F['ILR_i'][-1]:.4f} = {G['ILR_i'][-1]/F['ILR_i'][-1]:.3f}")
    print(f"  t=6  I3 qubit g/f = {G['I3_q'][-1]:.4f}/{F['I3_q'][-1]:.4f} | inv g/f = {G['I3_i'][-1]:.4f}/{F['I3_i'][-1]:.4f}")
    print("  I_inv(L:R)(t) gauge", np.round(G['ILR_i'], 3).tolist())
    print("  I_inv(L:R)(t) livre", np.round(F['ILR_i'], 3).tolist(), flush=True)
