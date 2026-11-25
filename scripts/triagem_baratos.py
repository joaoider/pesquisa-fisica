"""Triagem de observáveis baratos (v0.32): gauge × controles locais V∈{1,2,4} de mesmo gap, N=16.
O1 ν; O2 ⟨Q_iQ_j⟩_c (dist. 6, centro); O3 S₂ bloco central de 4; O4 F = Var(Q da meia cadeia) (= Var(E_meio) pela lei de Gauss).
Exato (média t∈[2,6]) + protocolo (Trotter dt=0,25, desfase p=0,01, K disparos na base computacional)."""
import json, time
import numpy as np
from pesquisa_fisica.emergence.ed_quench import evolve, gap, local_interaction_hamiltonian, gap_of, evolve_H, rdm
from pesquisa_fisica.emergence.protocolo_t3 import diag_energies, trotter_trajectory

N = 16
idx = np.arange(1 << N, dtype=np.int64)
bits = ((idx[:, None] >> np.arange(N)) & 1)
Qs = bits - (np.arange(N) % 2)                        # carga por sítio
Qleft = Qs[:, : N // 2].sum(1)
exc = np.where(np.arange(N) % 2 == 0, bits, 1 - bits).sum(1) / N
i0, j0 = 5, 11                                         # distância 6 em torno do centro
blk = [6, 7, 8, 9]


def diag_obs(p):
    """Observáveis diagonais a partir da distribuição p(config)."""
    nu = (p * exc).sum()
    qq = (p * Qs[:, i0] * Qs[:, j0]).sum() - (p * Qs[:, i0]).sum() * (p * Qs[:, j0]).sum()
    F = (p * Qleft ** 2).sum() - (p * Qleft).sum() ** 2
    return nu, qq, F


def S2(full):
    r = rdm(full, N, blk)
    return float(-np.log(np.real(np.trace(r @ r))))


def match(V, target):
    lo, hi = 0.0, 8.0
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if gap_of(local_interaction_hamiltonian(N, 1.0, mid, V)[0]) < target else (lo, mid)
    return 0.5 * (lo + hi)


times = np.arange(0, 6.0001, 0.5)
w = times >= 2
rec_steps = list(range(8, 25, 2))                      # t = 2.0, 2.5, ..., 6.0 com dt=0.25
out = []
rng = np.random.default_rng(7)
for (mg, J) in [(0.0, 1.0), (0.5, 2.0)]:
    g = gap(N, mg, J)
    models = {"gauge": ("g", mg, J)}
    for V in (1.0, 2.0, 4.0):
        models[f"V={V:g}"] = ("c", match(V, g), V)
    case = {"m": mg, "J": J, "gap": g, "models": {}}
    for name, (kind, mm, par) in models.items():
        t0 = time.time()
        if kind == "g":
            st = evolve(N, mm, par, times); E = diag_energies(N, mm, J=par)
        else:
            H, b = local_interaction_hamiltonian(N, 1.0, mm, par); st = evolve_H(H, b, N, times); E = diag_energies(N, mm, V=par)
        o = np.array([diag_obs(np.abs(s) ** 2) for s in st])
        s2 = np.array([S2(s) for s in st])
        ex = {"nu": float(o[w, 0].mean()), "QQ6": float(o[w, 1].mean()), "F": float(o[w, 2].mean()), "S2_blk4": float(s2[w].mean())}
        # protocolo: Trotter + desfase, 200 trajetórias → distribuição p em cada tempo registrado
        R = 200
        P = {k: np.zeros(1 << N) for k in rec_steps}
        for _ in range(R):
            tr = trotter_trajectory(N, E, 0.25, 24, 0.01, rng, record=set(rec_steps))
            for k in rec_steps:
                P[k] += np.abs(tr[k]) ** 2 / R
        Fprot = [diag_obs(P[k])[2] for k in rec_steps]
        # disparos: K=1e4 por tempo, 200 repetições do estimador da média temporal de F
        K, reps = 10000, 200
        est = []
        for _ in range(reps):
            vals = []
            for k in rec_steps:
                pk = P[k] / P[k].sum()
                smp = rng.choice(1 << N, size=K, p=pk)
                q = Qleft[smp]
                vals.append(q.var())
            est.append(np.mean(vals))
        case["models"][name] = {"m": mm, "par": par, "exato": ex, "F_protocolo": float(np.mean(Fprot)),
                                "F_est_media": float(np.mean(est)), "F_est_sigma": float(np.std(est)),
                                "secs": round(time.time() - t0, 1)}
        print(f"({mg},{J}) {name:6s} m={mm:.3f} | exato: ν={ex['nu']:.4f} QQ6={ex['QQ6']:+.5f} S2={ex['S2_blk4']:.3f} "
              f"F={ex['F']:.4f} | protocolo F={np.mean(Fprot):.4f} ± {np.std(est):.4f} ({case['models'][name]['secs']}s)", flush=True)
    out.append(case)
    json.dump(out, open("results/triagem_baratos_N16.json", "w"))
