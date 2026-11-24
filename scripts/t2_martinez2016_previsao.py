"""T2 (validação): previsão IDEAL para Martinez et al., Nature 534, 516 (2016) — N=4, J=w=1, m=0,5, ε₀=0, vácuo nu.
ν(t) = fração de sítios excitados (= (1/2N) Σ⟨(−1)^l σᶻ_l + 1⟩ na convenção do artigo);
E_N(t) = negatividade logarítmica (log₂) entre metades {1,2}|{3,4}. Calculado ANTES de ler as figuras."""
import json
import numpy as np
from pesquisa_fisica.emergence.ed_quench import evolve, rdm

N, m, J = 4, 0.5, 1.0
times = np.round(np.linspace(0, 3 * np.pi / 2, 121), 5)
st = evolve(N, m, J, times)


def nu(full):
    p = np.abs(full) ** 2
    idx = np.arange(1 << N)
    bits = (idx[:, None] >> np.arange(N)) & 1
    exc = np.where(np.arange(N) % 2 == 0, bits, 1 - bits)       # par ocupado = partícula; ímpar vazio = antipartícula
    return float((p[:, None] * exc).sum() / N)


def log_neg(full):
    # estado puro: E_N = 2 log₂ Σ √λ (λ = espectro de Schmidt entre metades)
    lam = np.clip(np.linalg.eigvalsh(rdm(full, N, [0, 1])), 0, None)
    return float(2 * np.log2(np.sqrt(lam).sum()))


out = {"params": {"N": N, "m_over_w": m, "J_over_w": J, "eps0": 0.0},
       "t": times.tolist(), "nu": [nu(s) for s in st], "E_N": [log_neg(s) for s in st]}
json.dump(out, open("results/t2_martinez2016_previsao.json", "w"))
nu_ = np.array(out["nu"]); en = np.array(out["E_N"])
for t in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0):
    k = int(np.argmin(np.abs(times - t)))
    print(f"wt={times[k]:.2f}  ν={nu_[k]:.4f}  E_N={en[k]:.4f}")
print("máx ν =", round(nu_.max(), 4), "em wt =", round(times[nu_.argmax()], 3), "| máx E_N =", round(en.max(), 4), "em wt =", round(times[en.argmax()], 3))
