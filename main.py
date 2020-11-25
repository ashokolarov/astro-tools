from nbody import *
from constants import *
from tools import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('ggplot')

# Bodies to propagate
bodies = [EARTH, MOON]

# Initial conditions
vel_moon = vorbital(EARTH["mass"], EARTH_MOON_DISTANCE)
init = [  EARTH_MOON_DISTANCE,       # x-position
                            0,       # y-position
                            0,       # z-position
                            0,       # x-velocity
                     vel_moon,       # y-velocity
                            0        # z-velocity
]

N  = 40000 # Number of nodes
dt = 200. # Timestep [s]

sim = TwoBodyPropagator(init, bodies, N, dt)
sim.propagate()
sim.plot()
