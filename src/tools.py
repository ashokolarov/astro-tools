import numpy as np
from constants import *


def v_circ(m, R):
    """
    Compute orbital velocity needed for a circular orbit
        m - Mass of orbited body [kg]
        R - Distance to center of body [m]

    Return: orbital velocity [m/s]
    """
    return np.sqrt(G * m / R)

def get_e(ra, rp):
    """
    Calculate eccentricity
        ra - apoapsis distance [m]
        rp - periapsis distance [m]

    Return: eccentricity [-]
    """
    return (ra - rp) / (ra + rp)

def get_a(ra, rp):
    """
    Calculate semi-major axis
        ra - apoapsis distance [m]
        rp - periapsis distance [m]

    Return: semi-major axis [m]

    """
    return (rp + ra) / 2

def get_rp(a, e):
    """
    Calculate periapsis
        a - semi-major axis [m]
        e - eccentricity [-]
    Return: periapsis distance [m]
    """
    return a * (1 - e)

def get_ra(a, e):
    """
    Calculate apoapsis
        a - semi-major axis [m]
        e - eccentricity [-]
    Return: apoapsis distance [m]
        """
    return a * (1 + e)

def KOE_TO_CSV(m, a, e, i, w, om, M0, t0, t):
    """
    Compile Cartesian State Vector from Keplerian Orbit Elements
        m - Mass of orbited body [kg]
        a - Semi-major axis [m]
        e - Eccentricity, 0<=e<=1 [-]
        i - Inclination [deg]
        w - Argument of periapsis [deg]
        om - Longitude of ascending node [deg]
        M0 - Mean anomaly at epoch t0 [deg]
        t0 - Epoch time [s]
        t - Considered epoch [s]
    Return: np.array holding x, y, z, xdot, ydot, zdot
    """
    mu = G * m
    w *= r2d
    om *= r2d
    i *= r2d
    M0 *= r2d

    if t0 != t:
        dt = t - t0
        M = M0 + dt * np.sqrt(mu / a ** 3)
        if M > 2 * np.pi:
            M %= 2 * np.pi
    else:
        M = M0

    tolerance = 1e-5
    E = M
    # Newton-Rhapson method to calculate eccentric anomaly up to a certain tolerance
    while E - e * np.sin(E) - M > tolerance:
        E = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))

    arg1 = np.sqrt(1 + e) * np.sin(E / 2)
    arg2 = np.sqrt(1 - e) * np.cos(E / 2)
    v = 2 * np.arctan2(arg1, arg2)

    rc = a * (1 - e * np.cos(E))
    o = np.array([np.cos(v), np.sin(v), 0], dtype=np.float64) * rc
    odot = np.array([-np.sin(E), np.sqrt(1 - e * e) * np.cos(E), 0], dtype=np.float64) * np.sqrt(mu * a) / rc

    ox, oy, _ = o
    ox_dot, oy_dot, _ = odot

    rx = ox * (np.cos(w) * np.cos(om) - np.sin(w) * np.cos(i) * np.sin(om)) - oy * (
            np.sin(w) * np.cos(om) + np.cos(w) * np.cos(i) * np.sin(om))
    ry = ox * (np.cos(w) * np.sin(om) + np.sin(om) * np.cos(i) * np.cos(om)) + oy * (
            np.cos(w) * np.cos(i) * np.cos(om) - np.sin(w) * np.sin(om))
    rz = ox * np.sin(w) * np.sin(i) + oy * np.cos(w) * np.sin(i)

    rx_dot = ox_dot * (np.cos(w) * np.cos(om) - np.sin(w) * np.cos(i) * np.sin(om)) - oy_dot * (
            np.sin(w) * np.cos(om) + np.cos(w) * np.cos(i) * np.sin(om))
    ry_dot = ox_dot * (np.cos(w) * np.sin(om) + np.sin(om) * np.cos(i) * np.cos(om)) + oy_dot * (
            np.cos(w) * np.cos(i) * np.cos(om) - np.sin(w) * np.sin(om))
    rz_dot = ox_dot * np.sin(w) * np.sin(i) + oy_dot * np.cos(w) * np.sin(i)

    return np.array([rx, ry, rz, rx_dot, ry_dot, rz_dot], dtype=np.float64)
