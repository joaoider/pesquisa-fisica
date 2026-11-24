"""v0.19: IM entre blocos de tamanho ℓ_b após o quench (exato, N=20) — o cone balístico reaparece em escala de méson?"""
import json, time
import numpy as np
from pesquisa_fisica.emergence.ed_quench import evolve, gap, block_mi, symmetric_blocks, peak
from pesquisa_fisica.emergence.schwinger_dmrg import free_gap

N = 20
times = np.round(np.arange(0, 6.0001, 0.1), 3)
LB, RS = [1, 2, 3, 4], [2, 4, 6, 8]


def match(N, target, lo=0.0, hi=6.0):
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if free_gap(N, mid) < target else (lo, mid)
    return 0.5 * (lo + hi)


def run(m, J):
    st = evolve(N, m, J, times)
    res = {}
    for lb in LB:
        for r in RS:
            A, B = symmetric_blocks(N, lb, r)
            I = [block_mi(s, N, A, B) for s in st]
            ok, tp, Im = peak(times, I)
            res[f"{lb}_{r}"] = {"I": I, "peak": ok, "t_peak": tp, "I_max": Im}
    return res


out = []
for m, J in [(0.0, 1.0), (0.5, 2.0)]:
    t0 = time.time()
    g = gap(N, m, J)
    mf = match(N, g)
    row = {"N": N, "m": m, "J": J, "gap": g, "m_free": mf, "times": times.tolist(),
           "gauge": run(m, J), "free": run(mf, 0.0)}
    row["secs"] = round(time.time() - t0, 1)
    out.append(row)
    json.dump(out, open("results/meson_coarse_N20.json", "w"))
    print(f"m={m} J={J} gap={g:.3f} mf={mf:.3f} ({row['secs']}s)")
    for side in ("gauge", "free"):
        for lb in LB:
            print(f"  {side:5s} lb={lb}", [(r, row[side][f'{lb}_{r}']['peak'], round(row[side][f'{lb}_{r}']['t_peak'], 2),
                                          round(row[side][f'{lb}_{r}']['I_max'], 4)) for r in RS], flush=True)
