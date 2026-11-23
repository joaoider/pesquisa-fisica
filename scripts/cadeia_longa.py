"""Teste de cadeia longa (v0.14): razão assintótica das taxas de decaimento gauge/livre a mesmo gap.
Uso: python scripts/cadeia_longa.py [N1,N2,...] [J1,J2,...] [chi]"""
import json, sys, time, warnings, logging
import numpy as np
warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
from pesquisa_fisica.emergence.schwinger_dmrg import (gauge_profiles, free_profiles, match_free_mass, decay_rate_fit)

Ns = [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else "40,64").split(",")]
Js = [float(x) for x in (sys.argv[2] if len(sys.argv) > 2 else "0.25,0.5,1.0").split(",")]
chi = int(sys.argv[3]) if len(sys.argv) > 3 else 120
out = []
for N in Ns:
    i = N // 4
    js = list(range(i + 2, 3 * N // 4 + 1))
    ds = np.array(js) - i
    for J in Js:
        t = time.time()
        g = gauge_profiles(N, J, i, js, chi_max=chi)
        mf = match_free_mass(N, g["gap"])
        f = free_profiles(N, mf, i, js)
        row = {"N": N, "J": J, "gap": g["gap"], "m_free": mf, "gap_free": f["gap"]}
        for k in ["C", "D", "MI_f"]:
            row[f"k_{k}_gauge"] = decay_rate_fit(ds, g[k])
            row[f"k_{k}_free"] = decay_rate_fit(ds, f[k])
        row["k_MI_s_gauge"] = decay_rate_fit(ds, g["MI_s"])
        row["ratio_MI_f"] = row["k_MI_f_gauge"] / row["k_MI_f_free"]
        row["ratio_MI_s_vs_free"] = row["k_MI_s_gauge"] / row["k_MI_f_free"]
        row["P1_kC_over_kD_gauge"] = row["k_C_gauge"] / row["k_D_gauge"]
        row["P2_kD_over_kC_free"] = row["k_D_free"] / row["k_C_free"]
        row["secs"] = round(time.time() - t, 1)
        row["profiles"] = {"d": ds.tolist(), **{f"{k}_gauge": np.asarray(g[k]).tolist() for k in ["C", "D", "MI_f", "MI_s"]},
                           **{f"{k}_free": np.asarray(f[k]).tolist() for k in ["C", "D", "MI_f"]}}
        out.append(row)
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items() if k != "profiles"}, flush=True)
name = sys.argv[4] if len(sys.argv) > 4 else "cadeia_longa"
json.dump(out, open(f"results/{name}.json", "w"), indent=1)
