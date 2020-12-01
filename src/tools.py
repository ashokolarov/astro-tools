import numpy as np
from constants import * 

def vorbital(m, R):
    """
    Compute velocity needed to maintain a circular velocity at a distance
    from the center of the body.

    m - Mass of orbited body [kg]
    R - Distance to center of body [m]
    """
    
    return np.sqrt(G * m / R)

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
    """
    mu = G * m
    w *= r2d
    om *= r2d
    i *= r2d
    M0 *= r2d

    if t0 != t:
        dt = t - t0
        M = M0 + dt * np.sqrt(mu / a**3)
        if M > 2*np.pi:
            M %= 2*np.pi
    else:
        M = M0

    tolerance = 1e-5
    E = M
    # Newton-Rhapson method to calculate eccentric anomaly up to a certain tolerance
    while (E - e * np.sin(E) - M > tolerance):
        E = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))

    arg1 = np.sqrt(1 + e) * np.sin(E/2)
    arg2 = np.sqrt(1 - e) * np.cos(E/2)
    v = 2 * np.arctan2(arg1, arg2)

    rc = a * (1 - e * np.cos(E))
    o = np.array([np.cos(v), np.sin(v), 0], dtype=np.float64) * rc
    odot = np.array([-np.sin(E), np.sqrt(1 - e*e) * np.cos(E), 0], dtype=np.float64) * np.sqrt(mu*a) / rc

    ox, oy, _ = o
    ox_dot, oy_dot, _ = odot

    rx = ox*(np.cos(w)*np.cos(om) - np.sin(w)*np.cos(i)*np.sin(om)) - oy*(np.sin(w)*np.cos(om) + np.cos(w)*np.cos(i)*np.sin(om))
    ry = ox*(np.cos(w)*np.sin(om) + np.sin(om)*np.cos(i)*np.cos(om)) + oy*(np.cos(w)*np.cos(i)*np.cos(om) - np.sin(w)*np.sin(om))
    rz = ox*np.sin(w)*np.sin(i) + oy*np.cos(w)*np.sin(i)

    rx_dot = ox_dot*(np.cos(w)*np.cos(om) - np.sin(w)*np.cos(i)*np.sin(om)) - oy_dot*(np.sin(w)*np.cos(om) + np.cos(w)*np.cos(i)*np.sin(om))
    ry_dot = ox_dot*(np.cos(w)*np.sin(om) + np.sin(om)*np.cos(i)*np.cos(om)) + oy_dot*(np.cos(w)*np.cos(i)*np.cos(om) - np.sin(w)*np.sin(om))
    rz_dot = ox_dot*np.sin(w)*np.sin(i) + oy_dot*np.cos(w)*np.sin(i)

    return np.array([rx,ry,rz,rx_dot,ry_dot,rz_dot], dtype=np.float64)

def CSV_TO_KOE(m, u):
    mu = G * m
    
    r = u[:3]
    v = u[3:]

    h = np.cross(r,v)
    h_norm = np.linalg.norm(h)

    r_norm = np.linalg.norm(r)
    v_norm = np.linalg.norm(v)

    E = ((v_norm*v_norm)/2) - (mu/r)

    a = -mu / (2*E)
    e = (1 - ( (h*h) / (a * mu)))**0.5

    i = np.arccos(h[2]/h_norm)
    om = np.arctan2(h[0], -h[1])

    p = a * (1 - e * e)
    v = np.arctan2(np.sqrt(p/mu) * np.dot(v,r), p - r_norm)

    wv = np.arctan2(r[2]/np.sin(i), r[0]*np.cos(om) + r[1]*np.sin(om))
    w = wv - v

    EC_ANOMALY = 2*np.arctan((((1-e)/(1+e))**0.5) * np.tan(v/2))
    M = EC_ANOMALY - e * np.sin(EC_ANOMALY)

    return [a, e, i, w, om, M]












    

    

    
        
    
    
    
