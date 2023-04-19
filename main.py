import numpy as np
import matplotlib.pyplot as plt

from nbody import orbitsys, body
from planets import earth, moon
global G
G = 6.6743e-20


if __name__ == '__main__':

    sys_config = {
        'tspan': [0, 60*60*24*15],
        'N': 1000,
        'propagate': 0,
    }

    body_config = {}

    e = body(body_config, planet=1)

    bodies = [e, ]

    nbodysys = orbitsys(sys_config, bodies)

    # plt.figure(1)
    # ax = plt.axes(projection='3d')
    # ax.plot3D(sc.state[:, 0], sc.state[:, 1], sc.state[:, 2])
    # ax.plot3D(sc.m_state[:, 0], sc.m_state[:, 1], sc.m_state[:, 2])
    # # limits to see moon
    # plt.xlim([-moon['radius_earth'], moon['radius_earth']])
    # plt.ylim([-moon['radius_earth'], moon['radius_earth']])
    # # limits to see satellite
    # # plt.xlim([-earth['radius']+1000,earth['radius']+1000])
    # # plt.ylim([-earth['radius']+1000,earth['radius']+1000])

    # plt.figure(2)
    # plt.plot(sc.scm[:, 0], sc.scm[:, 1])
    # plt.grid()
    # plt.show()
