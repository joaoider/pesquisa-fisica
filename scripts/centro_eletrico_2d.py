"""v0.28: centro elétrico em 2+1D (Z₂, rede 4×4 aberta): S_ext = H_borda + S_dest para todos os retângulos."""
import json
import numpy as np
from pesquisa_fisica.emergence.z2_links import Z2Lattice, region_entropies

lat = Z2Lattice(4, 4)
lb = lat.link_bits(np.arange(1 << lat.P, dtype=np.int64))
rects = [(x0, x1, y0, y1) for x0 in range(4) for x1 in range(x0 + 1, 5) for y0 in range(4) for y1 in range(y0 + 1, 5)
         if not (x0 == 0 and x1 == 4 and y0 == 0 and y1 == 4)]
out = []
for h in (0.15, 0.33, 0.8):
    E, gap, psi = lat.ground_state(h)
    rows = []
    for r in rects:
        d = region_entropies(lat, psi, lat.rect_links(*r), lb)
        d["rect"] = r
        rows.append(d)
    n = np.array([d["n_cut"] for d in rows]); S = np.array([d["S_ext"] for d in rows]); He = np.array([d["H_edge"] for d in rows])
    fitS = np.polyfit(n, S, 1); fitH = np.polyfit(n, He, 1)
    r2 = lambda y, f: 1 - ((y - np.polyval(f, n)) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    frac = He / np.where(S > 1e-12, S, np.nan)
    summ = {"h": h, "E": E, "gap": gap, "slope_S": fitS[0], "b_S": fitS[1], "b_S_over_ln2": fitS[1] / np.log(2),
            "R2_S": r2(S, fitS), "slope_H": fitH[0], "R2_H": r2(He, fitH), "edge_frac_min": float(np.nanmin(frac)),
            "edge_frac_mean": float(np.nanmean(frac)), "S_dist_max": float(max(d["S_dist"] for d in rows)),
            "S_dist_min": float(min(d["S_dist"] for d in rows))}
    out.append({"summary": summ, "regions": rows})
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in summ.items()}, flush=True)
json.dump(out, open("results/centro_eletrico_2d.json", "w"))
