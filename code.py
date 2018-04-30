
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
        self.m = np.zeros(self.N)

        self.dm = 1.0


    def set_initial_conditions(self, x0, v0, m0):
        self.x[0] = x0
        self.v[0] = v0
        self.m[0] = m0


    def set_time(self, time):
        self.time = time
        self.N = int(self.time/self.dt)

        self.x = np.zeros(self.N)
        self.v = np.zeros(self.N)
        self.t = np.zeros(self.N)
        self.m = np.zeros(self.N)


    def set_dt(self, dt):
        self.dt = dt


    def set_dm(self, dm):
        self.dm = dm


    def set_diff_eq(self, f):
        self.diff_eq = f;


    def rk4_step(self, x0, v0, t0, mi):
        dt = self.dt

        a1 = self.diff_eq(x0, v0, t0, mi)
        v1 = v0

        x_half_1 = x0 + v1 * dt/2.0
        v_half_1 = v0 + a1 * dt/2.0

        a2 = self.diff_eq(x_half_1, v_half_1, t0+dt/2.0, mi)
        v2 = v_half_1

        x_half_2 = x0 + v2 * dt/2.0
        v_half_2 = v0 + a2 * dt/2.0

        a3 = self.diff_eq(x_half_2, v_half_2, t0+dt/2.0, mi)
        v3 = v_half_2

        x_end = x0 + v3 * dt
        v_end = v0 + a3 * dt

        a4 = self.diff_eq(x_end, v_end, t0 + dt, mi)
        v4 = v_end

        a_middle = 1.0/6.0 * (a1 + 2*a2 + 2*a3 + a4)
        v_middle = 1.0/6.0 * (v1 + 2*v2 + 2*v3 + v4)

        x_end = x0 + v_middle * dt
        v_end = v0 + a_middle * dt

        return x_end, v_end


    def solve(self):
        x = self.x; v = self.v; t = self.t; N = self.N; dt = self.dt; m = self.m; dm = self.dm
        for i in range(N-1):
            [x[i+1], v[i+1]] = self.rk4_step(x[i], v[i], t[i], m[i])
            t[i+1] = t[i] + dt
            m[i+1] = m[i] + dm

        return x, v, t, m


    def solve_2(self):
        x = self.x; v = self.v; t = self.t; N = self.N; dt = self.dt; m = self.m; dm = self.dm
        for i in range(N-1):
            [x[i+1], v[i+1]] = self.rk4_step(x[i], v[i], t[i], m[i])
            t[i+1] = t[i] + dt
            m[i+1] = m[i] + dm
            if( True ):
                m[i+1] = m[i] - dm

        return x, v, t, m


def show_plots(k, m, ps, steps, lw):
    """
    Just a wrapper functions for plotting the PHASE SPACE, ENERGIES and POSITION/VELOCITIES.
    """
    # Plot phase space
    plt.plot(v, x, 'r', linewidth=2)
    plt.plot([-ps,ps],[0,0],'--k')
    plt.plot([0,0],[-ps,ps],'--k')
    plt.title('Phase space')
    plt.xlabel('Velocity m/s')
    plt.ylabel('Position')
    plt.show()

    # 
    #k = 1 # [N/m]
    Ek = 0.5*m*(v**2)
    Ep = 0.5*k*(x**2)

    # 
    plt.subplot(2,1,1)
    plt.plot(t, Ek, 'r', linewidth=lw)
    plt.plot(t, Ep, 'b', linewidth=lw)
    plt.plot(t, Ek+Ep, '--g', linewidth=1)
    plt.plot([0,steps], [0,0], '--k')
    plt.legend(['Kinetic energy', 'Potential energy'])

    # 
    plt.subplot(2,1,2)
    plt.plot(t, x, 'g', linewidth=lw)
    plt.plot(t, v, 'y', linewidth=lw)
    plt.plot([0,steps], [0,0], '--k')
    plt.legend(['Position', 'Velocity'])

    plt.show()


# 
def diff_eq_1(x, v, t, m):
    """
    Differential equation for problem 1
    """
    m = 0.5
    k = 1.0

    a = -(1./.5)*x

    return a


# 
def diff_eq_2(x, v, t, m):
    """
    Differential equation for problem 2
    """
    m = 0.5     # []
    k = 1.0     # []
    b = 0.1     # []

    a = -(k/m)*x - (b/m)*v

    return a


# 
def diff_eq_4(x, v, t, m):
    """
    Differential equation for problem 4
    """
    m = 0.5     # []
    k = 1.0     # []
    F_D = 0.7   # N
    w_0 = np.sqrt(k/m)
    omega_D = 13.0/(8.0*w_0)

    a = (F_D*np.cos(omega_D*t)-k*x)/m

    return a


# 
def diff_eq_5(x, v, t, m):
    """
    Differential equation for problem 5
    """
    m = 0.5     # []
    k = 1.0     # []
    b = 0.1     # []
    F_D = 0.7   # N
    w_0 = np.sqrt(k/m)
    omega_D = 13.0/(8.0*w_0)

    a = (F_D*np.cos(omega_D*t)-k*x-b*v)/m

    return a


# 
def diff_eq_6(x, v, t, m):
    """
    Differential equation for problem 6
    """
    k = 0.475       # []
    b = 0.001       # []
    g = 9.81        # 

    a = -(b*v + k*x + g)/m

    return a


solver = Solver(20.0, 1e-2)

problem_to_solve = 6
if( problem_to_solve == 1 ):
    # Initialize and solve the equation introduced in problem 1
    solver.set_diff_eq(diff_eq_1)
    solver.set_initial_conditions(1.0, 0.0, 0.5)
    [x, v, t, m] = solver.solve();

    show_plots(1.0, m[0], 1.5, 20.0, 2)

elif( problem_to_solve == 2 ):
    # 
    solver.set_diff_eq(diff_eq_2)
    solver.set_initial_conditions(1.0, 0.0, 0.5)
    [x, v, t, m] = solver.solve()

    show_plots(1.0, m[0], 1.5, 20.0, 2)

elif( problem_to_solve == 4):
    # 
    solver.set_diff_eq(diff_eq_4)
    solver.set_time(200.0)
    solver.set_initial_conditions(2.0, 0.0, 0.5)
    [x, v, t, m] = solver.solve()

    show_plots(1.0, m[0], 3.0, 200.0, 1)

elif( problem_to_solve == 5):
    # 
    solver.set_diff_eq(diff_eq_5)
    solver.set_time(100.0)
    solver.set_initial_conditions(2.0, 0.0, 0.5)
    [x, v, t, m] = solver.solve();

    show_plots(1.0, m[0], 3.0, 100.0, 1)
 
elif( problem_to_solve == 6):
    solver.set_diff_eq(diff_eq_6)
    solver.set_time(3.0)
    solver.set_dt(1e-4)
    solver.set_dm(0.00055)
    solver.set_initial_conditions(0.001, 0.001, 0.00001)
    [x, v, t, m] = solver.solve();

    show_plots(1.0, m[0], 3.0, 100.0, 1)

else:
    print "Please input a valid problem numer: [1,2,4]"


