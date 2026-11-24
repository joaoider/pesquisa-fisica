"""v0.25: supressão de pares × alcance R do potencial linear truncado (gap casado), blocos de 4 sítios a distância D."""
import json, time
import numpy as np
from pesquisa_fisica.emergence.ed_quench import truncated_linear_hamiltonian, gap_of, evolve_H, symmetric_blocks
from pesquisa_fisica.emergence.gauge_algebra import algebra_mi

N, J, m_gauge = 20, 2.0, 0.5
times = np.round(np.arange(0, 6.0001, 0.25), 3)
w = times >= 3
RS = [2, 3, 4, 5, 6, 8, 10, 12, 19]
rs = [4, 6, 8, 10]
target = gap_of(truncated_linear_hamiltonian(N, 1.0, m_gauge, J, 19)[0])


def match_m(R, lo=0.0, hi=6.0, iters=28):
    if R == 19:
        return m_gauge
    g0 = gap_of(truncated_linear_hamiltonian(N, 1.0, 0.0, J, R)[0])
    if g0 >= target:
        return 0.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        g = gap_of(truncated_linear_hamiltonian(N, 1.0, mid, J, R)[0])
        lo, hi = (mid, hi) if g < target else (lo, mid)
    return 0.5 * (lo + hi)


out = []
for R in RS:
    t0 = time.time()
    m = match_m(R)
    H, basis = truncated_linear_hamiltonian(N, 1.0, m, J, R)
    g = gap_of(H)
    st = evolve_H(H, basis, N, times)
    row = {"R": R, "m": m, "gap": g, "gap_target": target, "times": times.tolist(), "I": {}, "Ibar": {}}
    for r in rs:
        A, B = symmetric_blocks(N, 4, r)
        I = [algebra_mi(s, N, A, B, "inv") for s in st]
        row["I"][str(r)] = I
        row["Ibar"][str(r)] = float(np.mean(np.array(I)[w]))
    row["secs"] = round(time.time() - t0, 1)
    out.append(row)
    json.dump(out, open("results/escala_corda_N20.json", "w"))
    print(f"R={R:2d} m={m:.3f} gap={g:.3f}/{target:.3f} Ī(r=4,6,8,10)=" +
          str([f"{row['Ibar'][str(r)]:.5f}" for r in rs]) + f" ({row['secs']}s)", flush=True)
