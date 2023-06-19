import numpy as np
import matplotlib.pyplot as plt
from planets import earth, moon
import sys
# from gekko import GEKKO

global pi
pi = 3.14159
global G
G = 6.6743e-20  # km^3/kg/s^2


def default_config():
    config = {
        'dt': 0,
        'tspan': [0, 100],
        'N': 0,
        'propagate': 1,
        'gen_opt': 0,
        'j2': 0,
        'aero': 0,
    }
    return config


def default_body(type):
    if type:
        body = { # default satellite config
                'state0': np.array([]),
                'coes0': np.array([]),  # [a,e,i,AOP,RAAN,TA]
                'propagate': 1,
                'mass': 500,  # kg
                'Cd': 2.2,
        }
    else:
        body = { # default planet config
            'name': 'earth',
            'mass': 5.97219e24,  # kg
            'mu': 0.3986004418e6,  # km^3/s^2
            'radius': 6378.14,  # km
            'j2': 1.08262668355e-3,
        }

    return body


class body:

    def __init__(self, config, sat=0, planet=0,cb=0):
        if sat == 1:
            self.config = default_body(1)
        elif planet == 1:
            self.config = default_body(0)

        self.cb = cb
        self.sat = sat
        self.planet = planet

        for key in config.keys():
            # overwrite defaults, add additional config
            self.config[key] = config[key]

        if not planet:
            if any(self.config['state0']):
                self.config['state_or_coes'] = 1
            else:
                self.config['state_or_coes'] = 0


class orbitsys:

    def __init__(self, config, bodies):
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
        else:
            print("Error - no step size or step count specified")
            sys.exit()


        self.bodies = bodies
        self.planets = []
        for body in bodies: # place cb and non-cb planets into external data structures
            if body.planet == 1 and body.cb == 1:
                self.cb = body
                self.bodies.remove(body)
            elif body.planet == 1 and body.cb == 0:
                self.planets.append(body)
                self.bodies.remove(body)

        for body in bodies: # place cb and non-cb planets into external data structures
            if body.planet == 1 and body.cb == 0:
                self.planets.append(body)
                self.bodies.remove(body)
        
        if self.config['propagate']:
            self.prop_planets()
            self.prop_sats()



    def prop_planets(self):
        for planet in self.planets:
            planet.state = np.zeros([self.N,len(planet.config['state0'])])
            planet.state[0,:] = planet.config['state0']

        for i in range(self.N-1):
            for planet in self.planets:
                planet.state[i+1,:] = self.rk45(self.planet_eoms,planet.state[i,:],i)

    def prop_sats(self):
        for body in self.bodies:
            body.state = np.zeros([self.N,len(body.config['state0'])])
            body.state[0,:] = body.config['state0']

        for i in range(self.N-1):
            for sat in self.bodies:
                sat.state[i+1, :] = self.rk45(self.sat_eoms, sat.state[i, :], i)

    def sat_eoms(self, state, i):
        # sc eom wrt to cb
        r = np.zeros([3])
        v = np.zeros([3])
        r = state[0:3]
        v = state[3:6]
        dxdt = np.zeros([6])
        dxdt[0:3] = v
        # print(np.linalg.norm(r-self.bodies[0].state[i,0:3]))
        cbgrav = - r * self.cb.config['mu'] / np.linalg.norm(r)**3

        
        obgrav = 0
        j2pert = np.zeros([3])

        if self.config['j2'] == 1:
            R = np.linalg.norm(r)
            j2pert[0] += -3/2 * self.cb.config['j2'] *  self.cb.config['mu'] * self.cb.config['radius']**2 * r[0] / R**5 * (1 - 5*r[2]**2/R**2)
            j2pert[1] += -3/2 * self.cb.config['j2'] *  self.cb.config['mu'] * self.cb.config['radius']**2 * r[1] / R**5 * (1 - 5*r[2]**2/R**2)
            j2pert[2] += -3/2 * self.cb.config['j2'] *  self.cb.config['mu'] * self.cb.config['radius']**2 * r[2] / R**5 * (3 - 5*r[2]**2/R**2)
        
        for planet in self.planets:
            r_co = state[0:3] - planet.state[i-1,0:3]
            obgrav += - planet.config['mu'] / np.linalg.norm(r_co)**3 * r_co
            if self.config['j2'] == 1:
                R_co = np.linalg.norm(r)
                j2pert[0] += -3/2 * planet.config['j2'] *  planet.config['mu'] * planet.config['radius']**2 * r_co[0] / R_co**5 * (1 - 5*r_co[2]**2/R_co**2)
                j2pert[1] += -3/2 * planet.config['j2'] *  planet.config['mu'] * planet.config['radius']**2 * r_co[1] / R_co**5 * (1 - 5*r_co[2]**2/R_co**2)
                j2pert[2] += -3/2 * planet.config['j2'] *  planet.config['mu'] * planet.config['radius']**2 * r_co[2] / R_co**5 * (3 - 5*r_co[2]**2/R_co**2)
                
            
        dxdt[3:6] = cbgrav + obgrav + j2pert
        return dxdt

    def planet_eoms(self, state, i):
        perts = 0
        for planet in self.planets:
            r_co = state[0:3] - planet.state[i,0:3] # shold equal zero when current body = planet, isn't for some reason, close enough
            # print(np.linalg.norm(r_co[0]))
            if np.linalg.norm(r_co) > planet.config['radius']: # only calculate perturbation if house is cool
                perts += - planet.config['mu'] / np.linalg.norm(r_co)**3 * r_co

        r_cb_planet = np.zeros([3])
        v_cb_planet = np.zeros([3])
        r_cb_planet = state[0:3]
        v_cb_planet = state[3:6]
        dxdt = np.zeros([6])
        dxdt[0:3] = v_cb_planet
        dxdt[3:6] = - r_cb_planet * self.cb.config['mu'] / np.linalg.norm(r_cb_planet)**3 + perts

        return dxdt

    def rk4(self, f, x, i):
        k1 = self.dt * f(x, i)
        k2 = self.dt * f(x+k1/2, i)
        k3 = self.dt * f(x+k2/2, i)
        k4 = self.dt * f(x+k3, i)
        return x + 1/6 * (k1 + 2*k2 + 2*k3 + k4)

    def rk45(self, f, x, i): # rk45 formula with fixed step size, implement adaptive later
        B = np.array([[0, 0, 0, 0, 0],
                     [2/9, 0, 0, 0, 0],
                     [1/12, 1/4, 0, 0, 0],
                     [69/128, -243/128, 135/64, 0, 0],
                     [-17/12, 27/4, -27/5, 16/15, 0],
                     [65/432, -5/16, 13/16, 4/27, 5/144]])
        CH = np.array([47/450, 0, 12/25, 32/225, 1/30, 6/25])
        k1 = self.dt * f(x, i)
        k2 = self.dt * f(x + B[1, 0]*k1, i)
        k3 = self.dt * f(x + B[2, 0]*k1 + B[2, 1]*k2, i)
        k4 = self.dt * f(x + B[3, 0]*k1 + B[3, 1]*k2 + B[3, 2]*k3, i)
        k5 = self.dt * f(x + B[4, 0]*k1 + B[4, 1]*k2 +
                         B[4, 2]*k3 + B[4, 3]*k4, i)
        k6 = self.dt * f(x + B[5, 0]*k1 + B[5, 1]*k2 +
                         B[5, 2]*k3 + B[5, 3]*k4 + B[5, 4]*k5, i)
        return x + CH[0]*k1 + CH[1]*k2 + CH[2]*k3 + CH[3]*k4 + CH[4]*k5 + CH[5]*k6

    # reference: https://stackoverflow.com/questions/75285363/infeasibilities-solution-not-found-gekko-error
    def optimize_trajectory(self):
        pass


