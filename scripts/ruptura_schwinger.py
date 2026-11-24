"""v0.26: ressonância de ruptura no Schwinger massivo (sem corte): Ī_inv(D) entre células de 2 sítios, quench exato N=20."""
import json, time
import numpy as np
from pesquisa_fisica.emergence.ed_quench import evolve, gap
from pesquisa_fisica.emergence.gauge_algebra import algebra_mi

N = 20
times = np.round(np.arange(0, 6.0001, 0.25), 3)
w = times >= 3
Ds = list(range(2, 15, 2))


def cells(D):
    k = D // 2
    i = (k - 1) // 2
    ip = (k - 1) - i
    return list(range(8 - 2 * i, 10 - 2 * i)), list(range(10 + 2 * ip, 12 + 2 * ip))


out = []
for m, J in [(0.5, 2.0), (0.5, 1.0), (1.0, 1.0), (1.0, 0.75)]:
    t0 = time.time()
    g = gap(N, m, J)
    st = evolve(N, m, J, times)
    I = {D: [algebra_mi(s, N, *cells(D), "inv") for s in st] for D in Ds}
    Ibar = {D: float(np.mean(np.array(I[D])[w])) for D in Ds}
    row = {"m": m, "J": J, "gap": g, "d_star": 2 * g / J, "times": times.tolist(),
           "I": {str(D): I[D] for D in Ds}, "Ibar": {str(D): Ibar[D] for D in Ds}, "secs": round(time.time() - t0, 1)}
    out.append(row)
    json.dump(out, open("results/ruptura_schwinger_N20.json", "w"))
    print(f"m={m} J={J} d*={2*g/J:.2f} Ī(D=2..14)=" + str([f"{Ibar[D]:.5f}" for D in Ds]) + f" ({row['secs']}s)", flush=True)
