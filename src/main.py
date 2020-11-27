from nbody import *
from constants import *
from tools import *

# Bodies to propagate
satellites = [EARTH, MOON, APOLLO]

# Apollo
alt = 30000e3
h = alt + EARTH['radius']
vel = vorbital(EARTH['mass'], h)

# Moon
vel_moon = vorbital(EARTH['mass'], EARTH_MOON_DISTANCE)

# Initial conditions
init = [  0, EARTH_MOON_DISTANCE, h,     # x-position
          0, 0, 0,     # y-position                  
          0, 0, 0,     # z-position
          0, 0, 0,     # x-velocity
          0, vel, vel_moon,     # y-velocity  
          0, 0, 0      # z-velocity  
]

T_TOTAL = 24 * 60 * 60  # [1 day]
dt = 5.                # Timestep [s]
N  = int(T_TOTAL / dt)  # Number of nodes

sim = NBodyPropagator(init, satellites, N, dt)
sim.propagate()
sim.plot()
