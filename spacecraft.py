import numpy as np
import matplotlib.pyplot as plt
from planets import earth, moon
import scipy as sp
from scipy import optimize as opt
import sys
from gekko import GEKKO

global pi
pi = 3.14159
global G
G = 6.6743e-20  # km^3/kg/s^2


def default_config():
    config = {
        'state0': np.array([]),
        'm_state0': np.array([]),
        'moon_prop': 1,
        'coes': [],
        'perts': [],
        'dt': 0,
        'tspan': [0, 100],
        'N': 0,
        'propagate': 1,
        'gen_opt': 0,
        'specs': {
            'Cd': 2.2,
            'area': (1e-3)**2/4,
            'mass': 0,
        }
    }
    return config


class spacecraft:

    def __init__(self, config):
        self.config = default_config()

        for key in config.keys():
            self.config[key] = config[key]

        self.dt = self.config['dt']
        self.N = self.config['N']
        if self.dt:
            self.N = (self.config['tspan'][1] -
                      self.config['tspan'][0])/self.dt + 1
            self.N = int(self.N)
        elif self.N:
            self.dt = (self.config['tspan'][1] -
                       self.config['tspan'][0])/self.N
            # print(self.dt)
        else:
            print("Error - no step size or step count specified")
            sys.exit()

        self.state = np.zeros([self.N, len(self.config['state0'])])
        self.m_state = np.zeros([self.N, len(self.config['m_state0'])])
        self.state[0, :] = self.config['state0']
        self.m_state[0, :] = self.config['m_state0']
        self.scm = np.zeros([self.N, len(self.config['state0'])])
        self.scm[0, :] = self.config['state0'] - self.config['m_state0']

        # if self.config['coes']:
        #     pass
        #     # calculate states
        # elif self.config['state0']:
        #     pass
        #     # calculate coes

        self.prop_moon()

        if self.config['propagate']:
            self.prop_sc()

        self.t = np.zeros([self.N, 1])
        self.t[:, 0] = np.linspace(
            self.config['tspan'][0], self.config['tspan'][1], self.N)
        if self.config['gen_opt']:
            self.control = np.zeros([self.N, 3])
            self.opt_traj = self.optimize_trajectory()

    def prop_moon(self):
        for i in range(self.N-1):
            self.m_state[i+1,
                         :] = self.rk45(self.m_eoms, self.m_state[i, :], i)

    def prop_sc(self):
        for i in range(self.N-1):
            self.state[i+1, :] = self.rk45(self.sc_eoms, self.state[i, :], i)
            self.scm[i+1, :] = self.state[i+1, :] - self.m_state[i+1, :]

    def sc_eoms(self, state, i):
        r_EM = np.zeros([3])
        r_EM = self.m_state[i, 0:3]
        # sc eom wrt to earth
        r = np.zeros([3])
        v = np.zeros([3])
        r = state[0:3]
        v = state[3:6]
        dxdt = np.zeros([6])
        dxdt[0:3] = v

        egrav = - r * earth['mu'] / np.linalg.norm(r)**3
        mgrav = - (r - r_EM) * moon['mu'] / np.linalg.norm(r - r_EM)**3
        dxdt[3:6] = egrav + mgrav
        return dxdt

    def m_eoms(self, m_state, i):
        r_EM = np.zeros([3])
        v_EM = np.zeros([3])
        r_EM = m_state[0:3]
        v_EM = m_state[3:6]
        dxmdt = np.zeros([6])
        dxmdt[0:3] = v_EM
        dxmdt[3:6] = - r_EM * earth['mu'] / np.linalg.norm(r_EM)**3
        return dxmdt

    def rk4(self, f, x, i):
        k1 = self.dt * f(x, i)
        k2 = self.dt * f(x+k1/2, i)
        k3 = self.dt * f(x+k2/2, i)
        k4 = self.dt * f(x+k3, i)
        return x + 1/6 * (k1 + 2*k2 + 2*k3 + k4)

    def rk45(self, f, x, i):
        B = np.array([[0, 0, 0, 0, 0],
                     [2/9, 0, 0, 0, 0],
                     [1/12, 1/4, 0, 0, 0],
                     [69/128, -243/128, 135/64, 0, 0],
                     [-17/12, 27/4, -27/5, 16/15, 0],
                     [65/432, -5/16, 13/16, 4/27, 5/144]])
        CH = np.array([47/450, 0, 12/25, 32/225, 1/30, 6/25])
        k1 = self.dt * f(x, i)
        k2 = self.dt * f(x+B[1, 0]*k1, i)
        k3 = self.dt * f(x + B[2, 0]*k1 + B[2, 1]*k2, i)
        k4 = self.dt * f(x + B[3, 0]*k1 + B[3, 1]*k2 + B[3, 2]*k3, i)
        k5 = self.dt * f(x + B[4, 0]*k1 + B[4, 1]*k2 +
                         B[4, 2]*k3 + B[4, 3]*k4, i)
        k6 = self.dt * f(x + B[5, 0]*k1 + B[5, 1]*k2 +
                         B[5, 2]*k3 + B[5, 3]*k4 + B[5, 4]*k5, i)
        return x + CH[0]*k1 + CH[1]*k2 + CH[2]*k3 + CH[3]*k4 + CH[4]*k5 + CH[5]*k6

    def syseoms(self, x, u, t):
        dxdt = np.zeros([12])
        dxdt[0:6] = x[0:6]

        egrav = - x[0:3] * earth['mu'] / np.linalg.norm(x[0:3])**3
        mgrav = - (x[0:3] - x[3:6]) * moon['mu'] / \
            np.linalg.norm(x[0:3] - x[3:6])**3
        dxdt[6:9] = egrav + mgrav + u
        dxdt[9:12] = - x[3:6] * earth['mu'] / np.linalg.norm(x[3:6])**3
        return dxdt

    def obj_fun(self, x, u, t):
        J = 0
        for k in range(self.N):
            h = t[k+1] - t[k]
            J = J + (h / 2) * ((t[k]**2 + np.linalg.norm(x[k, :])**2 + np.linalg.norm(u[k, :])**2) + (
                t[k+1]**2 + np.linalg.norm(x[k+1, :])**2 + np.linalg.norm(u[k+1, :])**2))
        return J

    # def optimize_trajectory(self):
    #     # initial guesses
    #     init_state = np.concatenate((self.state,self.m_state),axis=1)
    #     init_control = np.array(self.control)
    #     init_time = np.array(self.t)
    #     print(np.shape(init_state),np.shape(init_control),np.shape(init_time))

    #     xinit = np.concatenate((init_state,init_control,init_time),axis=1)
    #     print(np.shape([xinit]))

    #     cons = (
    #         {'type': 'eq', 'fun': lambda x:  x[k+1, 0:3] - x[k, 0:3] - (t[k+1] - t[k])/2 * (
    #             self.syseoms(x[k, :], u[k, :], t[k]) + self.syseoms(x[k+1, :], u[k+1, :], t[k+1]))},
    #         # {'type': 'eq', 'fun': lambda x: (np.linalg.norm(
    #         #     x[self.N, 0:3] - x[self.N, 3:6]) - self.config['final_moon_radius'])},  # g1
    #         # # g2
    #         # {'type': 'eq', 'fun': lambda x: np.dot(x[self.N, 0:3]-x[self.N, 3:6], x[self.N, 6:9]-x[self.N,9:12])},
    #         # {'type': 'eq', 'fun': lambda x: x[0, :] - np.concatenate(self.config['state0'], self.config['m_state0'])},  # x[0] = x0
    #         # {'type': 'ineq', 'fun': lambda x, u, t, k: 1},
    #         # {'type': 'ineq', 'fun': lambda x, u, t, k: 1},
    #         # {'type': 'ineq', 'fun': lambda x, u, t, k: 1},
    #         # {'type': 'ineq', 'fun': lambda x, u, t, k: 1}
    #     )
    #     res = opt.minimize(lambda x: self.obj_fun(x), xinit, method='SLSQP', constraints=cons)
    #     return 0

    def optimize_trajectory(self):
        m = GEKKO(remote=False)

        m.time = np.linspace(
            self.config['tspan'][0], self.config['tspan'][1], self.N)
        tf = m.FV(value=self.config['tspan'][1],
                  lb=self.config['tspan'][0], ub=self.config['tspan'][1]*2)

        # state variables
        xs = m.Var(value=self.config['state0'][0])
        ys = m.Var(value=self.config['state0'][1])
        zs = m.Var(value=self.config['state0'][2])
        vxs = m.Var(value=self.config['state0'][3])
        vys = m.Var(value=self.config['state0'][4])
        vzs = m.Var(value=self.config['state0'][5])
        xm = m.Var(value=self.config['m_state0'][0])
        ym = m.Var(value=self.config['m_state0'][1])
        zm = m.Var(value=self.config['m_state0'][2])
        vxm = m.Var(value=self.config['m_state0'][3])
        vym = m.Var(value=self.config['m_state0'][4])
        vzm = m.Var(value=self.config['m_state0'][5])
        # control variables
        ux = m.MV(value=0)
        uy = m.MV(value=0)
        uz = m.MV(value=0)
        # intermediate variables
        xms = m.Intermediate(xs-xm)
        yms = m.Intermediate(ys-ym)
        zms = m.Intermediate(zs-zm)
        vxms = m.Intermediate(vxs-vxm)
        vyms = m.Intermediate(vys-vym)
        vzms = m.Intermediate(vzs-vzm)
        Rs = m.Intermediate((xs**2+ys**2+zs**2)**0.5)
        Rm = m.Intermediate((xm**2+ym**2+zm**2)**0.5)
        Rms = m.Intermediate((xms**2+yms**2+zms**2)**0.5)
        umag = m.Intermediate((ux**2+uy**2+uz**2)**0.5)
        axs = m.Intermediate(
            ux - earth['mu'] / Rs**3 * xs - moon['mu'] / Rms**3 * xms)
        ays = m.Intermediate(
            uy - earth['mu'] / Rs**3 * ys - moon['mu'] / Rms**3 * yms)
        azs = m.Intermediate(
            uz - earth['mu'] / Rs**3 * zs - moon['mu'] / Rms**3 * zms)
        axm = m.Intermediate(- earth['mu'] / Rm**3 * xm)
        aym = m.Intermediate(- earth['mu'] / Rm**3 * ym)
        azm = m.Intermediate(- earth['mu'] / Rm**3 * zm)

        # EOMS
        m.Equations((xs.dt() == vxs, ys.dt() == vys, zs.dt() == vzs,))
        m.Equations((vxs.dt() == axs, vys.dt() == ays, vzs.dt() == azs))
        m.Equations((xm.dt() == vxm, ym.dt() == vym, zm.dt() == vzm,))
        m.Equations((vxm.dt() == axm, vym.dt() == aym, vzm.dt() == azm))

        final = np.zeros_like(m.time)
        final[-1] = 1.0
        final = m.Param(value=final)

        intial = np.zeros_like(m.time)
        intial[0] = 1.0
        intial = m.Param(value=intial)

        # Inequality Constraints
        m.Equation(Rs > earth['radius'])
        m.Equation(Rms > moon['radius'])
        m.Equation(umag >= 0)
        m.Equation(umag <= 5)

        # Equality Constraints
        m.Equation(Rms * final == self.config['final_moon_radius'])
        m.Equation((xms*vxms*final)+(yms*vyms*final)+(yms*vyms*final) == 0)

        # Objective Function

        m.Obj(umag**2 + Rms**2 + tf*final)
        m.options.IMODE = 6
        m.options.MAX_ITER = 500
        m.options.NODES = 4

        m.solve(disp=True)
