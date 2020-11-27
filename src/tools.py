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

def KOE_TO_CSV(m, a, e, w, om, i, M0, t0, t):
    """
    Compile Cartesian State Vector from Keplerian Orbit Elements
    m - Mass of orbited body [kg]
    a - Semi-major axis [m]
    e - Eccentricity, 0<=e<=1 [-]
    w - Argument of periapsis [deg]
    om - Longitude of ascending node [deg]
    i - Inclination [deg]
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
    """
    Compile Cartesian State Vector from Keplerian Orbit Elements
    m - Mass of orbited body [kg]
    u - State vector [x,y,z,xdot,ydot,zdot]
    """
    
    mu = G * m
    
    r = u[:3]
    v = u[3:]

    h = np.cross(r,v)
    K = np.array([0,0,1])
    nhat = np.cross(K, h)

    norm_r = np.linalg.norm(r)
    norm_v = np.linalg.norm(v)

    e_vec = (((norm_v*norm_v - mu/norm_r) * r) - (np.dot(r,v)*v)) / mu
    e = np.linalg.norm(e_vec)

    E = norm_v*norm_v/2 - mu/norm_r
    a = -mu / (2 * E)
    p = a * (1 - e*e)

    norm_h = np.linalg.norm(h)
    norm_n = np.linalg.norm(nhat)
    
    i = np.arccos(h[2] / norm_h)
    om = np.arccos(nhat[0] / norm_n)
    if nhat[1] < 0:
        om = 2 * np.pi - om
        
    w = np.arccos(np.dot(nhat, e_vec) / (norm_n * e))
    if e_vec[1] < 0:
        w = 2 * np.pi - w

    M = np.arccos(np.dot(e_vec, r) / (e * norm_r))
    if np.dot(r,v) < 0:
        M = 2 * np.pi - M

    return [a, e, w, om, i, M]
    
        
    
    
    
