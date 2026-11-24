"""v0.23: I₃_inv e I(L:R) após o quench — gauge × controle com interação LOCAL (V Σ n n) de mesmo gap (exato N=20)."""
import json, time
import numpy as np
from pesquisa_fisica.emergence.ed_quench import local_interaction_hamiltonian, gap_of, evolve_H
from pesquisa_fisica.emergence.gauge_algebra import algebra_entropy, algebra_I3

N = 20
times = np.round(np.arange(0, 6.0001, 0.5), 3)
L, M, R = list(range(0, 6)), list(range(6, 14)), list(range(14, 20))
A, B, C = [7, 8], [9, 10], [11, 12]
ref = {(r["m"], r["J"]): r for r in json.load(open("results/quench_invariante_N20.json"))}


def match_m(V, target, lo=0.0, hi=6.0, iters=28):
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        g = gap_of(local_interaction_hamiltonian(N, 1.0, mid, V)[0])
        lo, hi = (mid, hi) if g < target else (lo, mid)
    m = 0.5 * (lo + hi)
    return m, gap_of(local_interaction_hamiltonian(N, 1.0, m, V)[0])


out = []
for key, r in ref.items():
    for V in (1.0, 2.0, 4.0):
        t0 = time.time()
        m, g = match_m(V, r["gap"])
        H, basis = local_interaction_hamiltonian(N, 1.0, m, V)
        st = evolve_H(H, basis, N, times)
        I3i = [algebra_I3(s, N, A, B, C, "inv") for s in st]
        I3q = [algebra_I3(s, N, A, B, C, "qubit") for s in st]
        ILR = [algebra_entropy(s, N, [L], (), "inv") + algebra_entropy(s, N, [R], (), "inv")
               - algebra_entropy(s, N, [L, R], (), "inv") for s in st]
        S_half = [algebra_entropy(s, N, [list(range(10))], (), "qubit") for s in st]
        row = {"gauge_case": list(key), "V": V, "m": m, "gap": g, "gap_target": r["gap"], "times": times.tolist(),
               "I3_i": I3i, "I3_q": I3q, "ILR_i": ILR, "S_half": S_half, "secs": round(time.time() - t0, 1)}
        out.append(row)
        json.dump(out, open("results/controle_interacao_N20.json", "w"))
        print(f"caso gauge {key} gap={r['gap']:.3f} | V={V} m={m:.3f} gap={g:.3f} | t=6: I3_inv={I3i[-1]:.4f} "
              f"(gauge {r['gauge']['I3_i'][-1]:.4f}, livre {r['free']['I3_i'][-1]:.4f}) I3_q={I3q[-1]:.4f} "
              f"I_inv(L:R)={ILR[-1]:.4f} (gauge {r['gauge']['ILR_i'][-1]:.4f}, livre {r['free']['ILR_i'][-1]:.4f}) "
              f"S_meia={S_half[-1]:.3f} (gauge {r['gauge']['S_half_q'][-1]:.3f}, livre {r['free']['S_half_q'][-1]:.3f}) ({row['secs']}s)", flush=True)
