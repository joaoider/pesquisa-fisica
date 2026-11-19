"""Cálculo pesado: SU(3) em 1+1D com N = 8 sítios (24 modos, 2 704 156 estados), diagonalização exata.

Roda no GitHub Actions (workflow "Cálculo pesado"). Grava results/su3_n8.json com:
energias, verificações (férmions livres, singleto), matriz de informação mútua por J,
cadeia recuperada (sim/não), razões de localidade e veredito das predições pré-registradas.

Uso local:  python scripts/su3_n8.py [--N 8] [--J 0 0.5 2 5]
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pesquisa_fisica.emergence.geometry_from_entanglement import knn_graph, spearman, upper  # noqa: E402
from pesquisa_fisica.emergence.sun_lattice import SUNChain, site_mutual_information  # noqa: E402


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=8)
    ap.add_argument("--Nc", type=int, default=3)
    ap.add_argument("--m", type=float, default=0.5)
    ap.add_argument("--J", type=float, nargs="+", default=[0.0, 0.5, 2.0, 5.0])
    ap.add_argument("--out", default="results/su3_n8.json")
    args = ap.parse_args()

    N = args.N
    D = np.abs(np.subtract.outer(np.arange(N), np.arange(N))).astype(float)
    cadeia = np.isclose(D, 1.0)
    mid = N // 2
    results = {"N": N, "Nc": args.Nc, "m": args.m, "w": 1.0, "maquina": platform.platform(),
               "por_J": {}}
    t_all = time.time()
    for J in args.J:
        t0 = time.time()
        ch = SUNChain(N, Nc=args.Nc, w=1.0, m=args.m, J=J)
        log(f"J={J}: dim={ch.dim}, construído em {time.time() - t0:.1f}s")
        E, v = ch.ground_state()
        log(f"J={J}: E={E:.12f} ({time.time() - t0:.1f}s)")
        S2 = ch.total_color_casimir(v)
        I = site_mutual_information(ch, v)
        A = knn_graph(I, 2, mutual=True)
        cadeia_ok = bool((A & cadeia).sum() == A.sum() == cadeia.sum())
        r = {
            "energia": E,
            "casimir_total": S2,
            "energia_fermions_livres": ch.free_energy() if J == 0 else None,
            "cadeia_recuperada": cadeia_ok,
            "spearman": spearman(upper(-I), upper(D)),
            "razao_vizinhos_dist2": float(I[mid - 1, mid] / I[mid - 2, mid]),
            "razao_vizinhos_dist3": float(I[mid - 1, mid] / I[mid - 2, mid + 1]),
            "informacao_mutua": I.tolist(),
            "segundos": time.time() - t0,
        }
        results["por_J"][str(J)] = r
        log(f"J={J}: S²={S2:.2e} cadeia={cadeia_ok} r(d=3)={r['razao_vizinhos_dist3']:.1f} "
            f"r(d=2)={r['razao_vizinhos_dist2']:.1f} ({r['segundos']:.0f}s)")
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(results, indent=2))     # salva parcial a cada J

    por = results["por_J"]
    Js = [str(J) for J in args.J]
    r3 = [por[j]["razao_vizinhos_dist3"] for j in Js]
    results["veredito"] = {
        "P1_fermions_livres": abs(por[Js[0]]["energia"] - por[Js[0]]["energia_fermions_livres"]) < 1e-8
        if args.J[0] == 0 else None,
        "P1_singleto": all(abs(por[j]["casimir_total"]) < 1e-7 for j in Js),
        "P2_cadeia_para_todo_J": all(por[j]["cadeia_recuperada"] for j in Js),
        "P3_localidade_cresce_com_J": all(a < b for a, b in zip(r3, r3[1:])),
    }
    results["segundos_total"] = time.time() - t_all
    Path(args.out).write_text(json.dumps(results, indent=2))
    log("VEREDITO:", json.dumps(results["veredito"]))
    if not all(v for v in results["veredito"].values() if v is not None):
        log("Atenção: alguma predição não se confirmou (isso é resultado, não erro do script).")


if __name__ == "__main__":
    main()
