"""v0.18: escala do emaranhamento após o quench — S(ℓ), I(L:R) e I₃; gauge vs livre de mesmo gap (TDVP, spin).
Uso: python scripts/escala_emaranhamento.py N tmax "m:J,..." chi nome"""
import json, sys, time, warnings, logging
import numpy as np
warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
from tenpy.algorithms.tdvp import TwoSiteTDVPEngine
from tenpy.networks.mps import MPS
from pesquisa_fisica.emergence.schwinger_dmrg import schwinger_mpo_model, ground_and_first_excited, match_free_mass
from pesquisa_fisica.emergence.quench import block_entropy, cut_entropy, tripartite_info

N = int(sys.argv[1]); tmax = float(sys.argv[2])
pts = [tuple(float(x) for x in p.split(":")) for p in sys.argv[3].split(",")]
chi = int(sys.argv[4]); name = sys.argv[5]
c = N // 2


def evolve(m, J, dt=0.05):
    M = schwinger_mpo_model(N, 1.0, m, J)
    psi = MPS.from_product_state(M.lat.mps_sites(), ["down" if k % 2 == 0 else "up" for k in range(N)], bc="finite")
    eng = TwoSiteTDVPEngine(psi, M, {"dt": dt, "N_steps": 1, "trunc_params": {"chi_max": chi, "svd_min": 1e-10}})
    for _ in range(int(round(tmax / dt))):
        eng.run()
    return psi


def measure(psi):
    rows = []
    for ell in range(1, 11):
        a = c - ell // 2
        Sb = block_entropy(psi, a, ell)
        rows.append({"ell": ell, "a": a, "S_block": Sb,
                     "I_LR": cut_entropy(psi, a) + cut_entropy(psi, a + ell) - Sb})
    return {"S_ell": rows, "S_cuts": [float(x) for x in psi.entanglement_entropy()],
            "I3_s1": tripartite_info(psi, c - 1, 1), "I3_s2": tripartite_info(psi, c - 3, 2)}


out = []
for m, J in pts:
    t0 = time.time()
    gap = ground_and_first_excited(N, 1.0, m, J, chi_max=80)["gap"]
    mf = match_free_mass(N, gap)
    row = {"N": N, "m": m, "J": J, "gap": float(gap), "m_free": float(mf), "tmax": tmax,
           "gauge": measure(evolve(m, J)), "free": measure(evolve(mf, 0.0))}
    row["secs"] = round(time.time() - t0, 1)
    out.append(row)
    for side in ("gauge", "free"):
        x = row[side]
        print(f"m={m} J={J} {side:5s} S(ℓ)={[round(r['S_block'],2) for r in x['S_ell']]} "
              f"I_LR={[round(r['I_LR'],3) for r in x['S_ell']]} I3_1={x['I3_s1']:.4f} I3_2={x['I3_s2']:.4f}", flush=True)
    json.dump(out, open(f"results/{name}.json", "w"))
