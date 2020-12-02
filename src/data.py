import numpy as np

EARTH_MOON_DISTANCE = 384400e3 # [m]

SUN = {
    "name" : "Sun",
    "radius" : 696340e3,   # [m]
    "mass" : 1.989e30      # [kg]
}

EARTH = {
    "name" : "Earth",
    "radius" : 6378e3,     # [m]
    "mass" : 5.972e24,     # [kg]
    "rot_vec" : np.array([0, 0, 7.2921e-5])  # [rad/s]
}

MOON = {
    "name" : "Moon",
    "radius" : 1737.1e3,   # [m]
    "mass" : 7.346e22      # [kg]
}

ISS = {
    "name" : "ISS",
    "radius" : 0,          # [m]
    "mass" : 419700,       # [kg]
    "Cd" : 3,              # [-]
    "A" : 200              # [m^2]
}

IRMAK = {
    "name" : "Irmak",
    "radius" : 0,          # [m]
    "mass" : 50 ,          # [kg]
    "Cd" : 3,              # [-]
    "A" : 200              # [m^2]
}

