from nbody import *
from data import *
from tools import *

# Bodies to propagate
satellites = [IRMAK]
koe = [EARTH['radius'] + 5000e3, 0.1, 100, 200, 90, 0, 0, 0]


# Initial conditions
init = KOE_TO_CSV(EARTH['mass'], *koe)


T_TOTAL = 4 * 60 * 60  
dt = 5.                # Timestep [s]
N  = int(T_TOTAL / dt)  # Number of nodes

sim = TwoBodyPropagator(init, EARTH, satellites, N, dt)
sim.propagate()
sim.plot()

