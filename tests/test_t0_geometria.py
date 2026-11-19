"""T0 — Consistência interna da ferramenta geométrica.

Antes de testar a física, testamos a régua: curvaturas conhecidas de métricas conhecidas.
"""
import sympy as sp

from pesquisa_fisica.geometry import Metric, flrw_flat, schwarzschild, sphere


def test_espaco_plano_tem_curvatura_nula():
    x, y, z = sp.symbols("x y z")
    m = Metric(sp.diag(1, 1, 1), (x, y, z))
    assert m.ricci_scalar == 0
    assert m.riemann_squared == 0


def test_esfera_2d_R_igual_2_sobre_r2():
    m = sphere(2)
    r = sp.Symbol("r", positive=True)
    assert sp.simplify(m.ricci_scalar - 2 / r**2) == 0


def test_esfera_3d_R_igual_6_sobre_r2():
    m = sphere(3)
    r = sp.Symbol("r", positive=True)
    assert sp.simplify(m.ricci_scalar - 6 / r**2) == 0


def test_schwarzschild_e_ricci_plana():
    m = schwarzschild()
    assert all(sp.simplify(e) == 0 for e in m.ricci)


def test_schwarzschild_kretschmann():
    m = schwarzschild()
    M, r = sp.Symbol("M", positive=True), sp.Symbol("r", real=True)
    assert sp.simplify(m.riemann_squared - 48 * M**2 / r**6) == 0


def test_flrw_escalar_de_ricci():
    t = sp.Symbol("t", real=True)
    a = sp.Function("a")(t)
    m = flrw_flat(a)
    esperado = 6 * (sp.diff(a, t, 2) / a + (sp.diff(a, t) / a) ** 2)
    assert sp.simplify(m.ricci_scalar - esperado) == 0
