import numpy as np

def create_sphere(R, coord):
    """
    Compute coordinates of a sphere centered at [0,0,0] using parametric equations.
        R - Radius of sphere [m]
        coord - X,Y,Z coordinates of the sphere's center
    Return:
        x,y,z - arrays holding coordinates of the generated sphere
    """
    x_offest,y_offset,z_offset = coord
    u = np.linspace(0, 2 * np.pi, 200)
    v = np.linspace(0, np.pi, 200)

    # The Cartesian coordinates of the unit sphere
    x = R * np.outer(np.cos(u), np.sin(v)) + x_offest
    y = R * np.outer(np.sin(u), np.sin(v)) + y_offset
    z = R * np.outer(np.ones(np.size(u)), np.cos(v)) + z_offset

    return x, y, z
