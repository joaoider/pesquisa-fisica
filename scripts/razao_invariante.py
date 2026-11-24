"""v0.21: razão gauge/livre (mesmo gap) da IM entre células, nas álgebras qubit, inv e carga (DMRG)."""
import json, sys, time, warnings, logging
import numpy as np
warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
from pesquisa_fisica.emergence.schwinger_dmrg import ground_and_first_excited, match_free_mass, decay_rate_power_fit
from pesquisa_fisica.emergence.quench import block_pair_rdms
from pesquisa_fisica.emergence.gauge_algebra import two_cell_mi

N = int(sys.argv[1]); Js = [float(x) for x in sys.argv[2].split(",")]; chi = int(sys.argv[3])
a0 = N // 4 if (N // 4) % 2 == 0 else N // 4 + 1
bs = list(range(a0 + 2, 3 * N // 4 + 1, 2))
MODES = ("qubit", "inv", "carga")


def profile(psi):
    rh = block_pair_rdms(psi, a0, bs, size=2)
    return {md: [two_cell_mi(rh[b], (a0, b))[md] for b in bs] for md in MODES}


out = []
for J in Js:
    t0 = time.time()
    g = ground_and_first_excited(N, 1.0, 0.0, J, chi_max=chi)
    mf = match_free_mass(N, g["gap"])
    f = ground_and_first_excited(N, 1.0, mf, 0.0, chi_max=chi, excited=False)
    pg, pf = profile(g["psi0"]), profile(f["psi0"])
    d = np.array(bs) - a0
    row = {"N": N, "J": J, "gap": float(g["gap"]), "m_free": float(mf), "d": d.tolist(), "gauge": pg, "free": pf}
    for md in MODES:
        kg, kf = decay_rate_power_fit(d, pg[md], floor=1e-12, parity=None), decay_rate_power_fit(d, pf[md], floor=1e-12, parity=None)
        row[f"k_{md}"] = {"gauge": kg, "free": kf, "ratio": kg / kf, "g_over_gap": kg / g["gap"], "f_over_gap": kf / g["gap"]}
    row["secs"] = round(time.time() - t0, 1)
    out.append(row)
    print(f"J={J} gap={g['gap']:.3f} " + " | ".join(
        f"{md}: κg={row[f'k_{md}']['gauge']:.3f} κf={row[f'k_{md}']['free']:.3f} razão={row[f'k_{md}']['ratio']:.3f} κg/gap={row[f'k_{md}']['g_over_gap']:.3f} κf/gap={row[f'k_{md}']['f_over_gap']:.3f}"
        for md in MODES) + f" ({row['secs']}s)", flush=True)
    json.dump(out, open(f"results/razao_invariante_N{N}.json", "w"))
