"""v0.24: potencial linear truncado em R (sem gauge), gap casado ao Schwinger; quench exato N=20, álgebra invariante."""
import json, time
import numpy as np
from pesquisa_fisica.emergence.ed_quench import truncated_linear_hamiltonian, gap_of, evolve_H
from pesquisa_fisica.emergence.gauge_algebra import algebra_entropy, algebra_I3

N = 20
times = np.round(np.arange(0, 6.0001, 0.5), 3)
L, R_, = list(range(0, 6)), list(range(14, 20))
A, B, C = [7, 8], [9, 10], [11, 12]
ref = {(r["m"], r["J"]): r for r in json.load(open("results/quench_invariante_N20.json"))}


def match_m(J, R, target, lo=0.0, hi=6.0, iters=28):
    g0 = gap_of(truncated_linear_hamiltonian(N, 1.0, 0.0, J, R)[0])
    if g0 >= target - 1e-9:
        return 0.0, g0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        g = gap_of(truncated_linear_hamiltonian(N, 1.0, mid, J, R)[0])
        lo, hi = (mid, hi) if g < target else (lo, mid)
    m = 0.5 * (lo + hi)
    return m, gap_of(truncated_linear_hamiltonian(N, 1.0, m, J, R)[0])


out = []
w = times >= 3
for (m0, J), r in ref.items():
    for R in (2, 3, 5, 8, 19):
        t0 = time.time()
        if R == 19:
            m, g = m0, r["gap"]
        else:
            m, g = match_m(J, R, r["gap"])
        H, basis = truncated_linear_hamiltonian(N, 1.0, m, J, R)
        st = evolve_H(H, basis, N, times)
        I3 = [algebra_I3(s, N, A, B, C, "inv") for s in st]
        ILR = [algebra_entropy(s, N, [L], (), "inv") + algebra_entropy(s, N, [R_], (), "inv")
               - algebra_entropy(s, N, [L, R_], (), "inv") for s in st]
        row = {"gauge_case": [m0, J], "R": R, "m": m, "gap": g, "gap_target": r["gap"], "times": times.tolist(),
               "I3_i": I3, "ILR_i": ILR, "I3_avg": float(np.mean(np.array(I3)[w])), "secs": round(time.time() - t0, 1)}
        out.append(row)
        json.dump(out, open("results/alcance_linear_N20.json", "w"))
        print(f"caso ({m0},{J}) R={R:2d} m={m:.3f} gap={g:.3f}/{r['gap']:.3f} | I3 médio={row['I3_avg']:.4f} "
              f"I3(t=6)={I3[-1]:.4f} I(L:R)(t=6)={ILR[-1]:.4f} ({row['secs']}s)", flush=True)
