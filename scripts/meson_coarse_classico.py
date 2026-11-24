"""v0.19 (diagnóstico): IM total × IM só de cargas (base σᶻ), blocos ℓ_b=4, exato N=20, t ≤ 6."""
import json
import numpy as np
from pesquisa_fisica.emergence.ed_quench import evolve, block_mi, classical_block_mi, symmetric_blocks, arrival_time

N = 20
times = np.round(np.arange(0, 6.0001, 0.1), 3)
R = json.load(open("results/meson_coarse_N20.json"))
out = []
for row in R:
    for side, (m, J) in (("gauge", (row["m"], row["J"])), ("free", (row["m_free"], 0.0))):
        st = evolve(N, m, J, times)
        rec = {"m": row["m"], "J": row["J"], "side": side, "r": {}}
        for r in (2, 4, 6, 8):
            A, B = symmetric_blocks(N, 4, r)
            Ic = [classical_block_mi(s, N, A, B) for s in st]
            rec["r"][str(r)] = {"I_class": Ic, "I_full": row[side][f"4_{r}"]["I"]}
        out.append(rec)
        ta = [arrival_time(times, rec["r"][str(r)]["I_class"], 1e-4) for r in (2, 4, 6, 8)]
        tf = [arrival_time(times, rec["r"][str(r)]["I_full"], 1e-3) for r in (2, 4, 6, 8)]
        print(row["m"], row["J"], side, "chegada total", np.round(tf, 2).tolist(), "| só cargas", np.round(ta, 2).tolist(), flush=True)
json.dump({"times": times.tolist(), "rows": out}, open("results/meson_coarse_classico_N20.json", "w"))
