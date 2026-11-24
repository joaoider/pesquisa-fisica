"""Regime confinante (v0.15): Schwinger m ≠ 0 vs livre a mesmo gap; taxa local κ_C(d) da correlação carregada.
Uso: python scripts/confinante.py N "m1:J1,m2:J2,..." chi nome"""
import json, sys, time, warnings, logging
import numpy as np
warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
from pesquisa_fisica.emergence.schwinger_dmrg import gauge_profiles, free_profiles, match_free_mass

N = int(sys.argv[1])
pts = [tuple(float(x) for x in p.split(":")) for p in sys.argv[2].split(",")]
chi = int(sys.argv[3]); name = sys.argv[4]
i = N // 4
js = list(range(i + 2, 3 * N // 4 + 1))
ds = np.array(js) - i
out = []
for m, J in pts:
    t = time.time()
    g = gauge_profiles(N, J, i, js, m=m, chi_max=chi)
    mf = match_free_mass(N, g["gap"])
    f = free_profiles(N, mf, i, js)
    row = {"N": N, "m": m, "J": J, "gap": g["gap"], "m_free": mf, "secs": None,
           "d": ds.tolist(), **{f"{k}_gauge": np.asarray(g[k]).tolist() for k in ["C", "D", "MI_f", "MI_s"]},
           **{f"{k}_free": np.asarray(f[k]).tolist() for k in ["C", "D", "MI_f"]}}
    row["secs"] = round(time.time() - t, 1)
    out.append(row)
    C = np.array(g["C"]); ev = (ds % 2 == 0) & (C > 1e-12)
    loc = -np.diff(np.log(C[ev])) / 2
    print(f"m={m} J={J} gap={g['gap']:.4f} mf={mf:.4f} t={row['secs']}s  kC_local(d={ds[ev][:-1][:12].tolist()})={np.round(loc[:12],3).tolist()}", flush=True)
    json.dump(out, open(f"results/{name}.json", "w"))
