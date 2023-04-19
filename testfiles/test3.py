from skyfield.api import load
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

ts = load.timescale()
tx = ts.utc(2010, 1, 1)
planets = load('de421.bsp')
earthx = planets['earth']
for keys in earthx:
    print(keys)