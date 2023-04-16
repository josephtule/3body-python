import numpy as np
import matplotlib.pyplot as plt

from spacecraft import spacecraft as SC
from planets import earth, moon
global G; G = 6.6743e-20


if __name__ == '__main__':
    sc_vc = np.sqrt(earth['mu']/ (earth['radius']+300))
    m_vc = np.sqrt(earth['mu']/ (moon['radius_earth']))
    config = {'state0': [0,-earth['radius']+300,0,sc_vc,0,0],
              'm_state0': [moon['radius_earth'],0,0,0,m_vc,0],
              'tspan': [0,60*60*24*31],
              'dt': 100
              }
    sc = SC(config)

    ax = plt.axes(projection='3d')
    ax.plot3D(sc.state[:,0],sc.state[:,1],sc.state[:,2])
    ax.plot3D(sc.m_state[:,0],sc.m_state[:,1],sc.m_state[:,2])
    ## limits to see moon
    # plt.xlim([-moon['radius_earth'],moon['radius_earth']])
    # plt.ylim([-moon['radius_earth'],moon['radius_earth']])
    # limits to see satellite
    plt.xlim([-earth['radius']+1000,earth['radius']+1000])
    plt.ylim([-earth['radius']+1000,earth['radius']+1000])
    plt.show()





     