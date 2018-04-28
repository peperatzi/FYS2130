
import numpy as np
import matplotlib.pyplot as plt


def diffEq(xi, vi, ti):
    #ai = f(xi, vi, ti)

    ai = f(xi, vi, ti)
    return ai

def rk(x0, v0, t0):
    a1 = diffEq(x1, v0, t0)
    v1 = v0

    x_half_1 = x0 + v1 * dt/2.0
    v_half_1 = v0 + a1 * dt/2.0

    a2 = diffEq(x_half_1, v_half_1, t0+dt/2.0)
    v2 = v_half_1

    x_half_2 = x0 + v2 * dt/2.0
    v_half_2 = v0 + a2 * dt/2.0

    a3 = diffEq(x_half_2, v_half_2, t0+dt/2.0
    v3 = v_half_2

    x_end = x0 + v3 * dt
    v_end = v0 + a3 * dt

    a4 = diffEq(x_end, v_end, t0 + dt)
    v4 = v_end

    a_middle = 1.0/6.0 * (a1 + 2*a2 + 2*a3 + a4)
    v_middle = 1.0/6.0 * (v1 + 2*v2 + 2*v3 + v4)

    x_end = x_start + v_middle * dt
    v_end = v_start + a_middle * dt

    return x_end, v_end



time = 10.0
dt = 1e-2

N = int(time/dt)

y = np.zeros(N)
V = np.zeros(N)
t = np.zeros(N)

for i in range(N):
    [y[i+1], v[i+1], t[i+1]] = rk(y[i], v[i], t[i])

