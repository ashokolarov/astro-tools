import numpy as np

EARTH_MOON_DISTANCE = 384400e3  # [m]

SUN = {
    "name": "Sun",
    "radius": 696340e3,  # [m]
    "mass": 1.989e30,  # [kg]
    "cmap": "autumn"
}

EARTH = {
    "name": "Earth",
    "radius": 6378e3,  # [m]
    "mass": 5.972e24,  # [kg]
    "rot_vec": np.array([0, 0, 7.2921e-5]),  # [rad/s]
    "J2": 1.08262668e-3,  # [-]
    "cmap": "GnBu"
}

MOON = {
    "name": "Moon",
    "radius": 1737.1e3,  # [m]
    "mass": 7.346e22,  # [kg]
    "cmap": "Greys"
}

ISS = {
    "name": "ISS",
    "mass": 419700,  # [kg]
    "Cd": 0.1,  # [-]
    "A": 20,  # [m^2]
    "Thrust" : 1.5e4 # [N]
}

GEO = {
    "name": "Geostationary satellite",
    "mass": 2000,  # [kg]
    "Cd": 1.5,  # [-]
    "A": 20,  # [m^2]
    "Thrust" : 500 # [N]
}

APOLLO = {
    "name": "Apollo",
    "mass": 28801, # [kg]
    "Cd": 0, # [-]
    "A": 0, # [m^2]
    "Thrust" : 91e3 # [N]
}

ENVI= {
    "name": "Envisat",
    "mass": 50,  # [kg]
    "Cd": 3,  # [-]
    "A": 2,  # [m^2]
    "Thrust" : 100 # [N]
}
