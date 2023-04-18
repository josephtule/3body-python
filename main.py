import numpy as np
import matplotlib.pyplot as plt

from spacecraft import spacecraft as SC
from planets import earth, moon
global G
G = 6.6743e-20


if __name__ == '__main__':
    moonalt = 500
    m_vc = np.sqrt(earth['mu'] / (moon['radius_earth']))
    sc_vc_moon = np.sqrt(moon['mu'] / (moon['radius']+moonalt)) + m_vc

    config_moonorbit = {
        'state0': np.array([moon['radius_earth']+moon['radius']+moonalt, 0, 0, 0, sc_vc_moon, 0]),
        'm_state0': np.array([moon['radius_earth'], 0, 0, 0, m_vc, 0]),
        'tspan': [0, 60*60*24*30],
        'dt': 100
    }

    earthalt = 250
    sc_v = np.sqrt(earth['mu'] * (2 / (earth['radius'] + earthalt) - 1/(earth['radius'] + earthalt)))
    config_hardlaunch = {
        'state0': np.array([0, -(earth['radius'] + earthalt), 0, sc_v*np.sqrt(2)*.99, 0, 0]),
        'm_state0': np.array([moon['radius_earth'], 0, 0, 0, m_vc, 0]),
        'tspan': [0, 60*60*24*30],
        'N': 10000,
        'propagate': 1,
    }

    plt.figure(1)
    sc = SC(config_hardlaunch)
    ax = plt.axes(projection='3d')
    ax.plot3D(sc.state[:, 0], sc.state[:, 1], sc.state[:, 2])
    ax.plot3D(sc.m_state[:, 0], sc.m_state[:, 1], sc.m_state[:, 2])
    # limits to see moon
    plt.xlim([-moon['radius_earth'], moon['radius_earth']])
    plt.ylim([-moon['radius_earth'], moon['radius_earth']])
    # limits to see satellite
    # plt.xlim([-earth['radius']+1000,earth['radius']+1000])
    # plt.ylim([-earth['radius']+1000,earth['radius']+1000])
    plt.show()

    plt.figure(2)
    plt.plot(sc.scm[:,0],sc.scm[:,1])
    plt.grid()
    plt.show()

    print(np.shape(sc.state))
