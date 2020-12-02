from nbody import *
from data import *

# Bodies to propagate
central_body = EARTH
satellites = [ISS]

# Keplerian orbital elements of sats
koe = [[EARTH['radius'] + 35000e3, 0.032, 87.87, 53.39, 227.9, 0.0]]

properties = {
       'J2' : True,
       'Drag' : False,
       'min_alt' : 200e3,
       'exit_upon_stopcondition' : False
}

T_TOTAL = 30 * 24 * 60 * 60
dt = 30.                # Time-step [s]
N = int(T_TOTAL / dt)   # Number of nodes

sim = TwoBodyPropagator.from_koe(koe, central_body, satellites, properties, N, dt)
sim.propagate()
sim.plot()
