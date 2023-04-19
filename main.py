import numpy as np
import matplotlib.pyplot as plt

from nbody import orbitsys, body
from planets import earth, moon

global G
G = 6.6743e-20

if __name__ == '__main__':

    sys_config = {
        'tspan': [0, 60*60*24*10],
        'N': 25000,
        'propagate': 1,
    }

    
    
    e = body(earth, planet=1,cb=1)
    m = moon; m['state0'] = [moon['radius_earth'],0,0,0,np.sqrt(earth['mu']/moon['radius_earth']),0]; 
    m = body(m,planet=1,cb=0)
    m2 = moon; m2['state0'] = [moon['radius_earth']*2,0,0,0,0,np.sqrt(earth['mu']/(moon['radius_earth']*2))]; m2['mu'] = 1e5
    m2 = body(m2,planet=1,cb=0)
    sc1 = body({'state0':[earth['radius'] + 250, 0, 0, 0, np.sqrt(earth['mu']/ (earth['radius'] + 250)), 0]},sat=1)
    sc2 = body({'state0':[earth['radius'] + 1000, 0, 0, 0, np.sqrt(earth['mu']/ (earth['radius'] + 1000)), 0]},sat=1)
    sc3 = body({'state0':[earth['radius'] + 1000, 0, 0, 0, 0, np.sqrt(earth['mu']/ (earth['radius'] + 1000))]},sat=1)
    bodies = [e, sc1, m]

    nbodysys = orbitsys(sys_config, bodies)
    
    plt.figure(1)
    ax = plt.axes(projection='3d')
    ax.plot3D(sc1.state[:, 0], sc1.state[:, 1], sc1.state[:, 2])
    # ax.plot3D(sc2.state[:, 0], sc2.state[:, 1], sc2.state[:, 2])
    # ax.plot3D(sc3.state[:, 0], sc3.state[:, 1], sc3.state[:, 2])
    ax.plot3D(m.state[:, 0], m.state[:, 1], m.state[:, 2])
    # ax.plot3D(m2.state[:, 0], m2.state[:, 1], m2.state[:, 2])
    # limits to see moon
    # ax.set_xlim([-moon['radius_earth'], moon['radius_earth']])
    # ax.set_ylim([-moon['radius_earth'], moon['radius_earth']])
    # ax.set_zlim([-moon['radius_earth'], moon['radius_earth']])

    # limits to see satellite
    ax.set_xlim([-earth['radius']+1000,earth['radius']+1000])
    ax.set_ylim([-earth['radius']+1000,earth['radius']+1000])
    ax.set_zlim([-earth['radius']+1000,earth['radius']+1000])

    plt.grid()
    plt.show()
