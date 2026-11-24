"""T2 — EXPLORATÓRIO (não pré-registrado): Trotter de 4 passos (dt=π/8) + desfase p por qubit por passo.
Compara com os pontos 'discretization errors' (vermelhos) e 'error model' (roxos) do próprio artigo e com os dados."""
import json
import numpy as np
from scipy.linalg import expm

N, w, m, J, p = 4, 1.0, 0.5, 1.0, 0.038
dim = 1 << N
idx = np.arange(dim)
bits = (idx[:, None] >> np.arange(N)) & 1
sz = 2 * bits - 1
stag = (-1.0) ** np.arange(N)
Q = (sz + stag) / 2
L = np.cumsum(Q, axis=1)[:, :-1]
Hz = np.diag(0.5 * m * (sz * stag).sum(1) + J * (L ** 2).sum(1))
def hop(bonds):
    H = np.zeros((dim, dim))
    for s in range(dim):
        for n in bonds:
            if ((s >> n) & 1) != ((s >> (n + 1)) & 1):
                H[s ^ ((1 << n) | (1 << (n + 1))), s] += w
    return H
Hx = hop(range(N - 1)); He, Ho = hop([0, 2]), hop([1])
H = Hz + Hx
exc = np.where(np.arange(N) % 2 == 0, bits, 1 - bits).sum(1) / N
psi0 = np.zeros(dim); psi0[sum(1 << k for k in range(1, N, 2))] = 1
rho0 = np.outer(psi0, psi0).astype(complex)
dt = np.pi / 8
Zs = [np.diag(np.where(bits[:, q] == 1, 1.0, -1.0)) for q in range(N)]


def dephase(r):
    for Z in Zs:
        r = (1 - p) * r + p * Z @ r @ Z
    return r


def run(U, noise):
    r, out = rho0.copy(), [0.0]
    for _ in range(4):
        r = U @ r @ U.conj().T
        if noise:
            r = dephase(r)
        out.append(float(np.real(np.trace(r @ np.diag(exc)))))
    return out


esquemas = {
    "ideal (exato)": expm(-1j * H * dt),
    "Trotter Z·X": expm(-1j * Hz * dt) @ expm(-1j * Hx * dt),
    "Trotter X·Z": expm(-1j * Hx * dt) @ expm(-1j * Hz * dt),
    "Trotter Z·par·ímpar": expm(-1j * Hz * dt) @ expm(-1j * Ho * dt) @ expm(-1j * He * dt),
}
D = json.load(open("results/t2_martinez2016_dados_digitalizados.json"))
red = [0.0, 0.1982, 0.4716, 0.4235, 0.2252]
print("artigo, discretização:", red)
print("artigo, modelo de erro:", D["nu_modelo_erro_artigo"])
print("dados:                 ", D["nu_exp"])
out = {}
for k, U in esquemas.items():
    a, b = run(U, False), run(U, True)
    out[k] = {"sem_ruido": a, "com_desfase": b}
    print(f"{k:22s} sem ruído {np.round(a,4).tolist()} | +desfase {np.round(b,4).tolist()} | "
          f"rms(dados) = {np.sqrt(np.mean((np.array(b)-D['nu_exp'])**2)):.4f}")
json.dump({"red_artigo": red, "esquemas": out}, open("results/t2_martinez2016_exploratorio.json", "w"), indent=1)
