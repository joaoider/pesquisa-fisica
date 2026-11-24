"""v0.20: IM em 4 álgebras (qubit, inv, ce, carga). Quench exato N=20 (m=0,J=1,ℓ_b=4) + vácuo N=16 (células)."""
import json
import numpy as np
from pesquisa_fisica.emergence.ed_quench import evolve, symmetric_blocks, peak, arrival_time
from pesquisa_fisica.emergence.gauge_algebra import algebra_mi
from pesquisa_fisica.emergence.schwinger import ground_state

MODES = ["qubit", "inv", "ce", "carga"]
out = {"quench": {}, "vacuo": {}}

N = 20
times = np.round(np.arange(0, 6.0001, 0.1), 3)
st = evolve(N, 0.0, 1.0, times)
out["quench"]["times"] = times.tolist()
for r in (2, 4, 6, 8):
    A, B = symmetric_blocks(N, 4, r)
    for md in MODES:
        I = [algebra_mi(s, N, A, B, md) for s in st]
        ok, tp, Im = peak(times, I)
        out["quench"][f"{md}_{r}"] = {"I": I, "peak": ok, "t_peak": tp, "I_max": Im}
for md in MODES:
    thr = 1e-4 if md in ("inv", "carga") else 1e-3
    ta = [arrival_time(times, out["quench"][f"{md}_{r}"]["I"], thr) for r in (2, 4, 6, 8)]
    v = 1 / (2 * np.polyfit([2, 4, 6, 8], ta, 1)[0])
    out["quench"][f"{md}_front"] = {"arrival": ta, "v": v, "thr": thr}
    print(md, "t_pico", [round(out['quench'][f'{md}_{r}']['t_peak'], 2) for r in (2, 4, 6, 8)],
          "I_max", [round(out['quench'][f'{md}_{r}']['I_max'], 4) for r in (2, 4, 6, 8)],
          "chegada", np.round(ta, 2).tolist(), f"v={v:.2f}", flush=True)

Nv = 16
for J in (0.0, 0.5, 2.0, 5.0):
    _, psi, _ = ground_state(Nv, w=1.0, m=0.5, J=J)
    A = [4, 5]
    row = {}
    for md in MODES:
        I1 = algebra_mi(psi, Nv, A, [6, 7], md)
        I3 = algebra_mi(psi, Nv, A, [10, 11], md)
        row[md] = {"I_d1": I1, "I_d3": I3, "ratio": I1 / I3 if I3 > 1e-15 else float("inf")}
    out["vacuo"][str(J)] = row
    print("vácuo J=", J, {md: (round(row[md]['I_d1'], 4), f"{row[md]['I_d3']:.2e}", round(row[md]['ratio'], 1)) for md in MODES}, flush=True)
json.dump(out, open("results/algebras_gauge.json", "w"))
