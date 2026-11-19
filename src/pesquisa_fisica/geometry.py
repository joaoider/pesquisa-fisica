"""Geometria diferencial simbólica: Christoffel, Riemann, Ricci e invariantes.

Convenções (MTW / Wald):
    Γ^a_{bc}   = 1/2 g^{ad} (∂_b g_{dc} + ∂_c g_{db} - ∂_d g_{bc})
    R^a_{bcd}  = ∂_c Γ^a_{db} - ∂_d Γ^a_{cb} + Γ^a_{ce} Γ^e_{db} - Γ^a_{de} Γ^e_{cb}
    R_{bd}     = R^a_{bad}
    R          = g^{bd} R_{bd}
Com essas convenções a esfera tem R > 0.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property

import sympy as sp


@dataclass
class Metric:
    g: sp.Matrix
    coords: tuple
    _simplify: bool = field(default=True, repr=False)

    def __post_init__(self):
        self.g = sp.Matrix(self.g)
        self.dim = len(self.coords)
        if self.g.shape != (self.dim, self.dim):
            raise ValueError("métrica e coordenadas com dimensões incompatíveis")

    def _s(self, expr):
        return sp.simplify(expr) if self._simplify else expr

    @cached_property
    def ginv(self) -> sp.Matrix:
        return self._s(self.g.inv())

    @cached_property
    def sqrt_abs_det(self):
        return sp.sqrt(sp.Abs(self._s(self.g.det())))

    @cached_property
    def christoffel(self):
        n, g, gi, x = self.dim, self.g, self.ginv, self.coords
        G = [[[0] * n for _ in range(n)] for _ in range(n)]
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    G[a][b][c] = self._s(sum(
                        gi[a, d] * (sp.diff(g[d, c], x[b]) + sp.diff(g[d, b], x[c]) - sp.diff(g[b, c], x[d]))
                        for d in range(n)) / 2)
        return G

    @cached_property
    def riemann(self):
        """R^a_{bcd}."""
        n, G, x = self.dim, self.christoffel, self.coords
        R = [[[[0] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    for d in range(n):
                        expr = sp.diff(G[a][d][b], x[c]) - sp.diff(G[a][c][b], x[d])
                        expr += sum(G[a][c][e] * G[e][d][b] - G[a][d][e] * G[e][c][b] for e in range(n))
                        R[a][b][c][d] = self._s(expr)
        return R

    @cached_property
    def riemann_lower(self):
        """R_{abcd} = g_{ae} R^e_{bcd}."""
        n, g, R = self.dim, self.g, self.riemann
        return [[[[self._s(sum(g[a, e] * R[e][b][c][d] for e in range(n)))
                   for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]

    @cached_property
    def ricci(self) -> sp.Matrix:
        n, R = self.dim, self.riemann
        return sp.Matrix(n, n, lambda b, d: self._s(sum(R[a][b][a][d] for a in range(n))))

    @cached_property
    def ricci_scalar(self):
        n = self.dim
        return self._s(sum(self.ginv[b, d] * self.ricci[b, d] for b in range(n) for d in range(n)))

    @cached_property
    def ricci_squared(self):
        """R_{ab} R^{ab}."""
        n, gi, Ric = self.dim, self.ginv, self.ricci
        up = gi * Ric * gi
        return self._s(sum(Ric[a, b] * up[a, b] for a in range(n) for b in range(n)))

    @cached_property
    def riemann_squared(self):
        """R_{abcd} R^{abcd} (escalar de Kretschmann)."""
        n, gi, Rl = self.dim, self.ginv, self.riemann_lower
        total = 0
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    for d in range(n):
                        if Rl[a][b][c][d] == 0:
                            continue
                        up = sum(gi[a, e] * gi[b, f] * gi[c, h] * gi[d, k] * Rl[e][f][h][k]
                                 for e in range(n) for f in range(n) for h in range(n) for k in range(n))
                        total += Rl[a][b][c][d] * up
        return self._s(total)

    def box(self, scalar):
        """□φ = (1/√|g|) ∂_a (√|g| g^{ab} ∂_b φ)."""
        n, x, gi, s = self.dim, self.coords, self.ginv, self.sqrt_abs_det
        return self._s(sum(sp.diff(s * gi[a, b] * sp.diff(scalar, x[b]), x[a]) for a in range(n) for b in range(n)) / s)


# ---------------------------------------------------------------- métricas de referência

def sphere(dim: int = 2, radius=None) -> Metric:
    """Esfera S^dim de raio r (dim = 2 ou 3)."""
    r = radius if radius is not None else sp.Symbol("r", positive=True)
    if dim == 2:
        th, ph = sp.symbols("theta phi", real=True)
        return Metric(sp.diag(r**2, r**2 * sp.sin(th)**2), (th, ph))
    if dim == 3:
        ch, th, ph = sp.symbols("chi theta phi", real=True)
        return Metric(sp.diag(r**2, r**2 * sp.sin(ch)**2, r**2 * sp.sin(ch)**2 * sp.sin(th)**2), (ch, th, ph))
    raise NotImplementedError


def schwarzschild(mass=None) -> Metric:
    M = mass if mass is not None else sp.Symbol("M", positive=True)
    t, r, th, ph = sp.symbols("t r theta phi", real=True)
    f = 1 - 2 * M / r
    return Metric(sp.diag(-f, 1 / f, r**2, r**2 * sp.sin(th)**2), (t, r, th, ph))


def flrw_flat(a=None) -> Metric:
    t, x, y, z = sp.symbols("t x y z", real=True)
    a = a if a is not None else sp.Function("a")(t)
    return Metric(sp.diag(-1, a**2, a**2, a**2), (t, x, y, z))
