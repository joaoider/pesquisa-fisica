"""v0.16: quench global (vácuo nu → H(m,J)); gauge (TDVP) vs livre exato de mesmo gap.
Uso: python scripts/quench_confinamento.py N tmax "m:J,..." chi nome"""
import json, sys, time, warnings, logging
import numpy as np
warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
from pesquisa_fisica.emergence.schwinger_dmrg import ground_and_first_excited, match_free_mass
from pesquisa_fisica.emergence.quench import gauge_quench, free_quench

N = int(sys.argv[1]); tmax = float(sys.argv[2])
pts = [tuple(float(x) for x in p.split(":")) for p in sys.argv[3].split(",")]
chi = int(sys.argv[4]); name = sys.argv[5]
times = np.round(np.arange(0, tmax + 1e-9, 0.5), 3)
c = N // 4
ds = list(range(1, N - c - 1))
out = []
for m, J in pts:
    t0 = time.time()
    gap = ground_and_first_excited(N, 1.0, m, J, chi_max=80)["gap"]
    mf = match_free_mass(N, gap)
    Sg, MIg, chis = gauge_quench(N, m, J, times, c, ds, chi_max=chi)
    Sf, MIf = free_quench(N, mf, times, c, ds)
    row = dict(N=N, m=m, J=J, gap=float(gap), m_free=float(mf), times=times.tolist(), ds=ds,
               S_gauge=Sg.tolist(), S_free=Sf.tolist(), MI_gauge=MIg.tolist(), MI_free=MIf.tolist(),
               chis=chis, secs=round(time.time() - t0, 1))
    out.append(row)
    print(f"m={m} J={J} gap={gap:.3f} mf={mf:.3f} S_g(tmax)={Sg[-1]:.3f} S_f(tmax)={Sf[-1]:.3f} "
          f"razao={Sg[-1]/Sf[-1]:.3f} chi={chis[-1]} t={row['secs']}s", flush=True)
    json.dump(out, open(f"results/{name}.json", "w"))
