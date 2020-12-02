import numpy as np
from constants import *


def KOE_TO_CSV(m, a, e, i, w, om, MA):
    """
    Compile Cartesian State Vector from Keplerian Orbit Elements
        m - Mass of orbited body [kg]
        a - Semi-major axis [m]
        e - Eccentricity, 0<=e<=1 [-]
        i - Inclination [deg]
        w - Argument of periapsis [deg]
        om - Longitude of ascending node [deg]
        MA - Mean Anomaly of the body [deg]
    Return: np.array holding x, y, z, xdot, ydot, zdot
    """
    mu = G * m
    i *= d2r
    w *= d2r
    om *= d2r

    n = np.sqrt(mu / a ** 3)

    tolerance = 1e-3
    EA = MA
    # Newton-Rhapson method to calculate eccentric anomaly up to a certain tolerance
    while EA - e * np.sin(EA) - MA > tolerance:
        EA = (EA + (MA - EA + e * np.sin(e)) / (1 - e * np.cos(e))) % (2 * np.pi)

    nu = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(EA / 2))

    r = a * (1 - e * np.cos(EA))
    h = np.sqrt(mu * a * (1 - e * e))

    X = r * (np.cos(om) * np.cos(w + nu) - np.sin(om) * np.sin(w + nu) * np.cos(i))
    Y = r * (np.sin(om) * np.cos(w + nu) + np.cos(om) * np.sin(w + nu) * np.cos(i))
    Z = r * np.sin(i) * np.sin(w + nu)

    p = a * (1 - e * e)

    Xdot = ((X * h * e) / (r * p)) * np.sin(nu) - (h / r) * (np.cos(om) * np.sin(w + nu) \
                                                            + np.sin(om) * np.cos(w + nu) * np.cos(i))

    Ydot = ((Y * h * e) / (r * p)) * np.sin(nu) - (h / r) * (np.sin(om) * np.sin(w + nu) \
                                                            - np.cos(om) * np.cos(w + nu) * np.cos(i))

    Zdot = ((Z * h * e) / (r * p)) * np.sin(nu) + (h / r) * np.sin(i) * np.cos(w + nu)

    return [X, Y, Z, Xdot, Ydot, Zdot]

def CSV_TO_KOE(m, u):
    mu = G * m
    r = u[:3]
    v = u[3:]

    r_norm = np.linalg.norm(r)
    v_norm = np.linalg.norm(v)

    h = np.cross(r, v)
    n = np.cross([0,0,1], h)
    e = np.cross(v, h)/mu - r/r_norm

    h_norm = np.linalg.norm(h)
    n_norm = np.linalg.norm(n)
    e_norm = np.linalg.norm(e)

    E = v_norm ** 2 / 2 - mu / r_norm
    a = -mu / (2 * E)
    i = np.arccos(h[2] / h_norm)
    om = np.arccos(n[0] / n_norm)
    if n[1] < 0:
        om = 2*np.pi - om
    w = np.arccos(np.dot(n, e) / (n_norm * e_norm))
    if e[2] < 0:
        w = 2*np.pi - w

    nu = np.arccos(np.dot(r,e) / (r_norm * e_norm))
    if np.dot(r,v) < 0:
        nu = 2*np.pi - nu

    EA = 2 * np.arctan(np.sqrt((1-e_norm)/(1+e_norm) * np.tan(nu/2)))
    MA = EA - e_norm * np.sin(EA)

    return [a,e_norm,i*r2d,w*r2d,om*r2d,MA]


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


"""
Coefficients of polynomials, used to determine air density at altitudes between 200km and 1000km.
http://www.braeunig.us/space/atmmodel.htm
"""
poly_coeff = [[1.1993e-9, -1.4511e-6, 6.9105e-4, -0.1736, -5.322],
              [1.1406e-10, -2.1308e-7, 1.5708e-4, -0.07029, -12.90],
              [8.1056e-12, -2.3584e-9, -2.6251e-6, -0.01563, -20.02],
              [-3.7012e-12, -8.6086e-9, 5.1188e-5, -0.06601, -6.138]]


def compute_poly(alt, lvl):
    """
    Compute polynomial of form:
    A × z^4 + B × z^3 + C × z^2 + D × z + E
    where coefficients A,B,C,D and E are to be found in poly_coeff array.
    """
    return poly_coeff[lvl][0] * alt ** 4 + poly_coeff[lvl][1] * alt ** 3 + poly_coeff[lvl][2] * alt ** 2 \
           + poly_coeff[lvl][3] * alt ** 1 + poly_coeff[0][4]


def earth_rho(alt):
    """
    Compute air density at at altitude of alt [m], according to http://www.braeunig.us/space/atmmodel.htm
        alt - altitude to calculate density at [m]

    Return:
        density [kg/m^3]
    """
    alt /= 1e3
    if 200 <= alt < 300:
        return np.exp(compute_poly(alt, 0))
    elif 300 <= alt < 500:
        return np.exp(compute_poly(alt, 1))
    elif 500 <= alt < 750:
        return np.exp(compute_poly(alt, 2))
    elif 750 <= alt < 1000:
        return np.exp(compute_poly(alt, 3))
    else:
        return 0
