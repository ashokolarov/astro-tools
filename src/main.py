from nbody import *
from data import *

# Bodies to propagate
central_body = EARTH
satellites = [ISS]

profile = None

# Keplerian orbital elements of sats
koe = [[EARTH['radius']+420e3, 0.0001866, 51.6414, 112.3890, 228.8997, 0]]

properties = {
       'J2' : True,
       'Drag' : True,
       'nbody': False,
       'min_alt' : 200e3,
       'exit_upon_stopcondition' : False
}

T_TOTAL = 30 * 24 * 60 * 60
dt = 20.                # Time-step [s]
sample_interval = 2 * 60 * 60

sim = TwoBodyPropagator.from_koe(koe, central_body, satellites, properties, profile, T_TOTAL, dt)
sim.propagate()
sim.plot_orbital_params(sample_interval)
