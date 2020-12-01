from nbody import *
from data import *
from tools import *

# Bodies to propagate
central_body = EARTH
satellites = [IRMAK, IRMAK, IRMAK]
# Keplerian orbital elements of sats
koe = [[EARTH['radius'] + 35000e3, 0.0, 0, 0, 0, 0, 0, 0],
       [EARTH['radius'] + 8000e3, 0.0, 20, 0, 0, 0, 0, 0],
       [EARTH['radius'] + 1000e3, 0.00, 0, 0, 0, 0, 0, 0]]

T_TOTAL = 24 * 60 * 60  
dt = 5.                # Timestep [s]
N  = int(T_TOTAL / dt)  # Number of nodes

sim = TwoBodyPropagator.from_koe(koe, central_body, satellites, N, dt)
sim.propagate()
sim.plot()
