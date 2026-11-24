"""v0.17: correlações multipartidas após o quench (N=24, t=6), gauge vs livre (ambos TDVP, partição de spin).
Uso: python scripts/multipartido.py N tmax "m:J,..." chi nome"""
import json, sys, time, warnings, logging
import numpy as np
warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
from tenpy.algorithms.tdvp import TwoSiteTDVPEngine
from tenpy.networks.mps import MPS
from pesquisa_fisica.emergence.schwinger_dmrg import schwinger_mpo_model, ground_and_first_excited, match_free_mass
from pesquisa_fisica.emergence.quench import block_pair_rdms, block_mi_and_excess, _ptrace, _S_rho

N = int(sys.argv[1]); tmax = float(sys.argv[2])
pts = [tuple(float(x) for x in p.split(":")) for p in sys.argv[3].split(",")]
chi = int(sys.argv[4]); name = sys.argv[5]


def evolve(m, J, dt=0.05):
    M = schwinger_mpo_model(N, 1.0, m, J)
    psi = MPS.from_product_state(M.lat.mps_sites(), ["down" if k % 2 == 0 else "up" for k in range(N)], bc="finite")
    eng = TwoSiteTDVPEngine(psi, M, {"dt": dt, "N_steps": 1, "trunc_params": {"chi_max": chi, "svd_min": 1e-10}})
    for _ in range(int(round(tmax / dt))):
        eng.run()
    return psi


def measure(psi):
    res = {}
    for label, a0 in (("alinhado", 4), ("desalinhado", 5)):
        bs = list(range(a0 + 2, N - 1, 2))
        rh = block_pair_rdms(psi, a0, bs, size=2)
        rows = []
        for b in bs:
            I_blk, I_pairs = block_mi_and_excess(rh[b])
            rows.append({"d": b - a0, "I_block": I_blk, "I_pairs": I_pairs})
        res[label] = rows
    a0 = 4
    rh1 = block_pair_rdms(psi, a0, list(range(a0 + 1, N)), size=1)
    res["site"] = [{"d": b - a0, "I": _S_rho(_ptrace(r, [0], 2)) + _S_rho(_ptrace(r, [1], 2)) - _S_rho(r)}
                   for b, r in rh1.items()]
    res["S_half"] = float(psi.entanglement_entropy()[N // 2 - 1])
    return res


out = []
for m, J in pts:
    t0 = time.time()
    gap = ground_and_first_excited(N, 1.0, m, J, chi_max=80)["gap"]
    mf = match_free_mass(N, gap)
    row = {"N": N, "m": m, "J": J, "gap": float(gap), "m_free": float(mf), "tmax": tmax,
           "gauge": measure(evolve(m, J)), "free": measure(evolve(mf, 0.0))}
    row["secs"] = round(time.time() - t0, 1)
    out.append(row)
    g, f = row["gauge"]["alinhado"], row["free"]["alinhado"]
    print(f"m={m} J={J} t={row['secs']}s", [(r['d'], round(r['I_block'], 4), round(r['I_pairs'], 4)) for r in g[:6]],
          "| livre", [(r['d'], round(r['I_block'], 4), round(r['I_pairs'], 4)) for r in f[:6]], flush=True)
    json.dump(out, open(f"results/{name}.json", "w"))
