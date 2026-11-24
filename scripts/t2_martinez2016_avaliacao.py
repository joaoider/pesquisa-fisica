"""T2 (validação): avalia P1–P3 pré-registrados contra os dados digitalizados de Martinez et al. 2016."""
import json
import numpy as np
from pesquisa_fisica.emergence.ed_quench import evolve, rdm

D = json.load(open("results/t2_martinez2016_dados_digitalizados.json"))
N, m, J = 4, 0.5, 1.0


def nu(full):
    p = np.abs(full) ** 2; idx = np.arange(1 << N); b = (idx[:, None] >> np.arange(N)) & 1
    return float((p[:, None] * np.where(np.arange(N) % 2 == 0, b, 1 - b)).sum() / N)


def EN(full):
    lam = np.clip(np.linalg.eigvalsh(rdm(full, N, [0, 1])), 0, None)
    return float(2 * np.log2(np.sqrt(lam).sum()))


# P1: curva ideal do artigo × nossa
pts = np.array(D["curva_ideal_artigo_px"])
t_art = (pts[:, 0] - 67.1) / (503.03 - 67.1) * np.pi / 2
nu_art = (285.24 - pts[:, 1]) / 554.4
nu_nos = np.array([nu(s) for s in evolve(N, m, J, t_art)])
dP1 = np.abs(nu_art - nu_nos)
# P2, P3: pontos experimentais
wt = np.array(D["wt"]); st = evolve(N, m, J, wt)
nu_id = np.array([nu(s) for s in st]); en_id = np.array([EN(s) for s in st])
dnu = np.array(D["nu_exp"]) - nu_id
en_exp = np.array(D["E_n_exp_m05"])
res = {"P1_max_dnu": float(dP1.max()), "P1": bool(dP1.max() <= 0.03),
       "nu_ideal": nu_id.tolist(), "dnu_exp_menos_ideal": dnu.tolist(), "P2": bool(np.all(np.abs(dnu) <= 0.07)),
       "E_N_ideal": en_id.tolist(), "E_n_exp": en_exp.tolist(),
       "P3": bool(np.all(en_exp[1:] > 0) and np.all(en_exp <= en_id + 0.1))}
json.dump(res, open("results/t2_martinez2016_avaliacao.json", "w"), indent=1)
print(f"P1: máx |ν_artigo − ν_nosso| ao longo da curva ideal = {dP1.max():.4f}  → {'PASSOU' if res['P1'] else 'FALHOU'}")
for k in range(5):
    print(f"  wt={wt[k]:.3f}: ν ideal {nu_id[k]:.4f} | exp {D['nu_exp'][k]:.4f} ± {D['nu_exp_err'][k]:.4f} | Δ={dnu[k]:+.4f}"
          f"   E_N ideal {en_id[k]:.3f} | exp {en_exp[k]:.3f}")
print("P2:", "PASSOU" if res["P2"] else "FALHOU", "| P3:", "PASSOU" if res["P3"] else "FALHOU")
