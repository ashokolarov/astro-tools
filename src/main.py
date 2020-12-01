from nbody import *
from data import *

# Bodies to propagate
central_body = EARTH
satellites = [MOON]

# Keplerian orbital elements of sats
koe = [[EARTH_MOON_DISTANCE, 0.0, 0, 0, 0, 0, 0, 0]]

pertubations = {
       'J2' : False,
}

T_TOTAL = 7 * 24 * 60 * 60
dt = 30.                # Timestep [s]
N  = int(T_TOTAL / dt)  # Number of nodes

sim = TwoBodyPropagator.from_koe(koe, central_body, satellites, pertubations, N, dt)
sim.propagate()
sim.plot()
