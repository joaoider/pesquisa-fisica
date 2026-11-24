"""v0.29: Z_N (N=2..5) em 2+1D, rede 3×3 aberta: fração destilável S_dest/S_ext × N (rumo a U(1))."""
import json, sys, time, gc
import numpy as np
from pesquisa_fisica.emergence.zn_links import ZNLattice, region_entropies_zn

h = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
Ns = [int(x) for x in (sys.argv[2] if len(sys.argv) > 2 else "2,3,4,5").split(",")]
rects = [(x0, x1, y0, y1) for x0 in range(3) for x1 in range(x0 + 1, 4) for y0 in range(3) for y1 in range(y0 + 1, 4)
         if not (x0 == 0 and x1 == 3 and y0 == 0 and y1 == 3)]
out = []
for N in Ns:
    t0 = time.time()
    lat = ZNLattice(3, 3, N)
    E0, gap, psi = lat.ground_state(h)
    E = lat.link_values(lat.digits()).astype(np.int64)
    rows = []
    for r in rects:
        d = region_entropies_zn(lat, psi, lat.rect_links(*r), E)
        if d is not None and d["S_ext"] > 1e-10:
            d["rect"] = r; rows.append(d)
    n = np.array([d["n_cut"] for d in rows]); He = np.array([d["H_edge"] for d in rows])
    frac = np.array([d["S_dist"] / d["S_ext"] for d in rows])
    fit = np.polyfit(n, He, 1)
    R2 = 1 - ((He - np.polyval(fit, n)) ** 2).sum() / ((He - He.mean()) ** 2).sum()
    summ = {"N": N, "h": h, "E0": E0, "gap": gap, "n_regions": len(rows), "frac_dist_mean": float(frac.mean()),
            "frac_dist_max": float(frac.max()), "R2_edge": float(R2), "slope_edge": float(fit[0]),
            "S_dist_mean": float(np.mean([d["S_dist"] for d in rows])), "secs": round(time.time() - t0, 1)}
    out.append({"summary": summ, "regions": rows})
    json.dump(out, open(f"results/zn_fracao_bulk_h{h}.json", "w"))
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in summ.items()}, flush=True)
    del lat, psi, E; gc.collect()
