#!/usr/bin/env python3
"""Independent check that the traveling-wave formula is an exact solution.

Evaluates the PDE residual

    R(x, t) = tau u_tt + u_t - D u_xx - M (u - u^3)

at random points for  u(x,t) = 1/(1+exp[lambda (x - c t)]),
lambda = sqrt(2+9 tau)/2, c = 3/sqrt(2+9 tau), tau = 0.2, D = M = 1,
using exact analytic derivatives. If max|R| is at the level of machine
precision, the formula is an exact (bounded) solution of the dimensional
PDE with the parameters used in the OOFEM input.
"""
import numpy as np

TAU, D, M = 0.2, 1.0, 1.0
LAMBDA = np.sqrt(2.0 + 9.0 * TAU) / 2.0
C = 1.5 / LAMBDA

rng = np.random.default_rng(2026)
x = rng.uniform(-12, 12, 2000)
t = rng.uniform(0, 0.2, 2000)
z = LAMBDA * (x - C * t)
u = 1.0 / (1.0 + np.exp(z))
# xi = x - c t;  u_xi = -lambda u(1-u);  u_xixi = lambda^2 u(1-u)(1-2u)
uxi = -LAMBDA * u * (1.0 - u)
uxxi = LAMBDA ** 2 * u * (1.0 - u) * (1.0 - 2.0 * u)
res = TAU * C ** 2 * uxxi - C * uxi - D * uxxi - M * (u - u ** 3)
print(f'lambda = {LAMBDA!r}')
print(f'speed c = {C!r}   (exact interface x_0(t) = c t)')
print(f'interface at T = 0.2: x_0 = {C * 0.2!r}')
print(f'max |PDE residual| over 2000 random points: {np.max(np.abs(res)):.3e}')
