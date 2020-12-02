from nbody import *
from data import *

# Bodies to propagate
central_body = EARTH
satellites = [ISS]

# Keplerian orbital elements of sats
koe = [[EARTH['radius'] + 3590e3, 0.032, 87.87, 53.39, 227.9, 0.0]]

properties = {
       'J2' : False,
       'Drag' : False,
       'min_alt' : 200e3,
       'exit_upon_stopcondition' : False
}

T_TOTAL = 24 * 60 * 60
dt = 5.                # Time-step [s]
N = int(T_TOTAL / dt)   # Number of nodes

sim = TwoBodyPropagator.from_koe(koe, central_body, satellites, properties, N, dt)
sim.propagate()
sim.plot()
