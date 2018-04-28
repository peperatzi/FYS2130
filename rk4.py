
import numpy as np
import matplotlib.pyplot as plt


class Solver():
    def __init__(self, time, dt):
        self.time = time
        self.dt = dt
        self.N = int(self.time/self.dt)

        self.x = np.zeros(self.N)
        self.v = np.zeros(self.N)
        self.t = np.zeros(self.N)


    def set_initial_conditions(self, x0, v0):
        self.x[0] = x0
        self.v[0] = v0


    def set_diff_eq(self, f):
        self.diff_eq = f;


    def rk4_step(self, x0, v0, t0):
        dt = self.dt

        a1 = self.diff_eq(x0, v0, t0)
        v1 = v0

        x_half_1 = x0 + v1 * dt/2.0
        v_half_1 = v0 + a1 * dt/2.0

        a2 = self.diff_eq(x_half_1, v_half_1, t0+dt/2.0)
        v2 = v_half_1

        x_half_2 = x0 + v2 * dt/2.0
        v_half_2 = v0 + a2 * dt/2.0

        a3 = self.diff_eq(x_half_2, v_half_2, t0+dt/2.0)
        v3 = v_half_2

        x_end = x0 + v3 * dt
        v_end = v0 + a3 * dt

        a4 = self.diff_eq(x_end, v_end, t0 + dt)
        v4 = v_end

        a_middle = 1.0/6.0 * (a1 + 2*a2 + 2*a3 + a4)
        v_middle = 1.0/6.0 * (v1 + 2*v2 + 2*v3 + v4)

        x_end = x0 + v_middle * dt
        v_end = v0 + a_middle * dt

        return x_end, v_end


    def solve(self):
        x = self.x; v = self.v; t = self.t; N = self.N; dt = self.dt
        for i in range(N-1):
            [x[i+1], v[i+1]] = self.rk4_step(x[i], v[i], t[i])
            t[i+1] = t[i] + dt 

        return x, v, t


# Initialze solver
solver = Solver(20.0, 1e-2)


# 
def diff_eq_1(x, v, t):
    """
    Differential equation for problem 1
    """
    m = 0.5
    k = 1.0
    a = -(1./.5)*x
    return a


# 
def diff_eq_2(x, v, t):
    """
    Differential equation for problem 2
    """
    m = 0.5     # []
    k = 1.0     # []
    b = 0.1     # []
    a = -(k/m)*x - (b/m)*v
    return a



# 
def diff_eq_4(x, v, t):
    """
    Differential equation for problem 2
    """
    m = 0.5     # []
    k = 1.0     # []
    b = 0.1     # []
    a = -(k/m)*x - (b/m)*v
    return a


# Prepare problem 1
solver.set_diff_eq(diff_eq_2)
solver.set_initial_conditions(1.0, 0.0)

# Yeah babe!
[x, v, t] = solver.solve();

# Plot this pure awsomness
plt.plot(x, v)
plt.show()



