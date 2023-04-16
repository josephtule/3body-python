import numpy as np
import matplotlib.pyplot as plt
from planets import earth, moon
import scipy as sp

global pi; pi = 3.14159
global G; G = 6.6743e-20  # km^3/kg/s^2



def default_config():
    config = {
        'state0': np.array([]),
        'm_state0': np.array([]),
        'moon_prop': 1,
        'coes': [],
        'perts': [],
        'dt': 1,
        'tspan': [0, 100],
        'N': 0,

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
        if self.dt:
            self.N = (self.config['tspan'][1] - self.config['tspan'][0])/self.dt + 1
            self.N = int(self.N)

   
        self.state = np.zeros([self.N,len(self.config['state0'])])
        self.m_state = np.zeros([self.N,len(self.config['m_state0'])])
        self.state[0,:] = self.config['state0']
        self.m_state[0,:] = self.config['m_state0']


        if self.config['coes']:
            pass
            # calculate states
        elif self.config['state0']:
            pass
            # calculate coes

        self.prop_moon()
        self.prop_sc()

    def prop_moon(self):
        for i in range(self.N-1):
            self.m_state[i+1,:] = self.rk4(self.m_eoms,self.m_state[i,:],i)

    def prop_sc(self):
        for i in range(self.N-1):
            self.state[i+1,:] = self.rk4(self.sc_eoms,self.state[i,:],i)

    def sc_eoms(self,state,i):
        r_EM = np.zeros([3])
        r_EM = self.m_state[i,0:3]
        # sc eom wrt to earth
        r = np.zeros([3]); v = np.zeros([3])
        r = state[0:3]; v = state[3:6]
        dxdt = np.zeros([6])
        dxdt[0:3] = v
        
        egrav = - r * earth['mu'] / np.linalg.norm(r)**3
        mgrav = - 0*(r - r_EM) * moon['mu'] / np.linalg.norm(r - r_EM)**3
        dxdt[3:6] = egrav + mgrav
        return dxdt

    def m_eoms(self,m_state,i):
        r_EM = np.zeros([3]); v_EM = np.zeros([3])
        r_EM = m_state[0:3]
        v_EM = m_state[3:6]
        dxmdt = np.zeros([6])
        dxmdt[0:3] = v_EM
        dxmdt[3:6] = - r_EM * earth['mu']/ np.linalg.norm(r_EM)**3
        return dxmdt

    def rk4(self,f,x,i):
        k1 = self.dt * f(x,i)
        k2 = self.dt * f(x+k1/2,i)
        k3 = self.dt * f(x+k2/2,i)
        k4 = self.dt * f(x+k3,i)
        return x + 1/6 * (k1 + 2*k2 + 2*k3 + k4)



    # def rk4(self,f, i):
    #     B = np.array([[0, 0, 0, 0, 0],
    #                  [2/9, 0, 0, 0, 0],
    #                  [1/12, 1/4, 0, 0, 0],
    #                  [69/128, -243/128, 135/64, 0, 0],
    #                  [-17/12, 27/4, -27/5, 16/15, 0],
    #                  [65/432, -5/16, 13/16, 4/27, 5/144]])
    #     CH = np.array([47/450, 0, 12/25, 32/225, 1/30, 6/25])
    #     k1 = self.dt * f(self.state(:, i))
    #     k2 = self.dt * f(self.state(: , i) + B(2, 1)*k1)
    #     k3 = self.dt * f(self.state(:, i) + B(3, 1)*k1 + B(3, 2)*k2)
    #     k4 = self.dt * f(self.state(: , i) + B(4, 1)*k1 + B(4, 2)*k2 + B(4, 3)*k3)
    #     k5 = self.dt * f(self.state(:, i) + B(5, 1)*k1 + B(5, 2)*k2 + B(5, 3)*k3 + B(5, 4)*k4)
    #     k6 = self.dt * f(self.state(: , i) + B(6, 1)*k1 + B(6, 2)*k2 + B(6, 3)*k3 + B(6, 4)*k4 + B(6, 5)*k5)
    #     self.state(:, i+1) = sc.state(: , i) + CH(1)*k1 + CH(2)*k2 + CH(3)*k3 + CH(4)*k4 + CH(5)*k5 + CH(6)*k6
