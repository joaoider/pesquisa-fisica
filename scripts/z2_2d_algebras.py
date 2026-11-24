"""v0.27: gauge Z₂ pura 2+1D (dual TFIM 5×4): IM entre plaquetes nas álgebras full × inv; geometria emergente."""
import json, time
import numpy as np
from pesquisa_fisica.emergence.z2_gauge_2d import ground_state, mi_matrices
from pesquisa_fisica.emergence.geometry_from_entanglement import (emergent_distances, spearman, true_distances, upper,
                                                                  classical_mds)

Lx, Ly = 5, 4
n = Lx * Ly
Dtrue = true_distances((Lx, Ly))
near = np.isclose(Dtrue, 1.0)
far = np.isclose(Dtrue, Dtrue.max())
out = []
for h in (0.15, 0.33, 0.8):
    t0 = time.time()
    E, gap, full = ground_state(Lx, Ly, h)
    If, Ii = mi_matrices(full, n)
    row = {"h": h, "E": E, "gap": gap, "I_full": If.tolist(), "I_inv": Ii.tolist()}
    for k, I in (("full", If), ("inv", Ii)):
        nr, fr = I[near].mean(), I[far].mean()
        D = emergent_distances(I)
        sp = spearman(upper(D), upper(Dtrue))
        ev = classical_mds(D)[0] if isinstance(classical_mds(D), tuple) else None
        row[k] = {"near": float(nr), "far": float(fr), "far_over_near": float(fr / nr), "spearman": sp}
    row["secs"] = round(time.time() - t0, 1)
    out.append(row)
    print(f"h={h} gap={gap:.3f} | full: perto={row['full']['near']:.4e} longe={row['full']['far']:.4e} "
          f"razão={row['full']['far_over_near']:.3e} Spearman={row['full']['spearman']:.3f} | inv: perto={row['inv']['near']:.4e} "
          f"longe={row['inv']['far']:.4e} razão={row['inv']['far_over_near']:.3e} Spearman={row['inv']['spearman']:.3f} ({row['secs']}s)", flush=True)
json.dump(out, open("results/z2_2d_algebras.json", "w"))
