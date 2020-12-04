import numpy as np
from constants import *
from tools import *
from plotter import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
mpl.rcParams['axes.formatter.useoffset'] = False



class Propagator():
    """
    Propagator class used to store general data about the simulation and propagate it.
    Params:
        n_body - number of bodies in the sim
        bodies - list of dictionary holding the bodies included in the sim
        N - number of timesteps
        dt - size of timestep
        T_tot - total time of simulation
    """

    def __init__(self, bodies, T_tot, dt):
        """
        Constructor
        Args:
            bodies - list of dictionaries, holds the bodies to be propagated
            T_tot - int, total runtime of simulation
            dt - float, size of timestep
        """
        self.n_body = len(bodies)
        self.bodies = bodies
        self.T_tot = T_tot
        self.dt = dt
        self.N = int(T_tot / dt)
        self.Tcur = 0
        self.run = True

    def rk4(self, f, u):
        """
        Runge-Kutta explicit 4-stage scheme - 1
          t  - time
          u  - solution at t
        Return: approximation of y at t+dt.
        """
        k1 = f(u)
        k2 = f(u + self.dt * k1 / 2)
        k3 = f(u + self.dt * k2 / 2)
        k4 = f(u + self.dt * k3)
        return u + self.dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

class TwoBodyPropagator(Propagator):

    def __init__(self, init, central_body, satellites, properties, profile, T_tot, dt):
        """
        Constructor
        """
        super().__init__(satellites, T_tot, dt)
        self.central = central_body
        self.properties = properties
        self.profile = profile
        self.M_Multibody = 6 * self.n_body
        self.u = np.zeros((self.N, self.M_Multibody), dtype=np.float64)
        self.u[0, :] = init

    @classmethod
    def from_koe(cls, koe, central_body, satellites, properties, profile, T_tot, dt):
        num_sat = len(satellites)
        states = []
        for i in range(num_sat):
            u = KOE_TO_CSV(central_body['mass'], *koe[i])
            states.append(u)
        init = np.array(states).flatten('F')
        return cls(init, central_body, satellites, properties, profile, T_tot, dt)

    def f_twobody(self, u):
        x = u[0:self.n_body]
        y = u[self.n_body:2 * self.n_body]
        z = u[2 * self.n_body:3 * self.n_body]
        xdot = u[3 * self.n_body:4 * self.n_body]
        ydot = u[4 * self.n_body:5 * self.n_body]
        zdot = u[5 * self.n_body:]
        f = np.zeros((self.n_body, 6))
        f[:, 0] = xdot
        f[:, 1] = ydot
        f[:, 2] = zdot
        for i in range(0, self.n_body):
            R = np.sqrt((x[i]) ** 2 + (y[i]) ** 2 + (z[i]) ** 2)
            alt = R - self.central['radius']
            f[i, 3] = -G * self.central['mass'] * x[i] / R ** 3
            f[i, 4] = -G * self.central['mass'] * y[i] / R ** 3
            f[i, 5] = -G * self.central['mass'] * z[i] / R ** 3

            if self.properties['J2'] == True:
                coeff = (3 * self.central['J2'] * self.central['mass'] * G * self.central['radius'] ** 2) / (2 * R ** 5)
                f[i, 3] += coeff * (5 * (z[i]**2 / R ** 2) - 1) * x[i]
                f[i, 4] += coeff * (5 * (z[i]**2 / R ** 2) - 1) * y[i]
                f[i, 5] += coeff * (5 * (z[i]**2 / R ** 2) - 3) * z[i]

            if self.properties['Drag'] == True:
                rho = earth_rho(alt)
                r = np.array([x[i], y[i], z[i]])
                v = np.array([xdot[i], ydot[i], zdot[i]])
                v_rel = v - np.cross(self.central['rot_vec'], r)
                speed = np.linalg.norm(v_rel)
                a_drag = -v_rel*0.5*rho*speed*self.bodies[i]['Cd']*self.bodies[i]['A'] / self.bodies[i]['mass']
                f[i, 3:6] += a_drag

            if self.properties['nbody'] == True:
                for j in range(0, self.n_body):
                    if i != j:
                        r = np.sqrt((x[j] - x[i]) ** 2 + (y[j] - y[i]) ** 2 + (z[j] - z[i]) ** 2)
                        f[i, 3] += G * self.bodies[j]['mass'] * (x[j] - x[i]) / r ** 3
                        f[i, 4] += G * self.bodies[j]['mass'] * (y[j] - y[i]) / r ** 3
                        f[i, 5] += G * self.bodies[j]['mass'] * (z[j] - z[i]) / r ** 3

            if self.profile is not None:
                if self.Tcur > self.profile['t0'] and self.Tcur < self.profile['t1'] and 'Thrust' in self.bodies[i]:
                    self.dt = self.profile['dt']
                    v = np.array([xdot[i], ydot[i], zdot[i]])
                    vhat = v / np.linalg.norm(v)
                    a_t = self.bodies[i]['Thrust'] * vhat / self.bodies[i]['mass']
                    f[i, 3:6] += a_t
                else:
                    self.dt = 10

            if alt < self.properties['min_alt'] and self.properties['exit_upon_stopcondition']:
                print(f"{self.bodies[i]['name']} has re-entered {self.central['name']} and crashed.")
                self.run = False

        return f.T.flatten()

    def propagate(self):
        for i in range(self.N - 1):
            if self.run:
                self.Tcur += self.dt
                self.u[i + 1, :] = self.rk4(self.f_twobody, self.u[i])
            else:
                self.u = self.u[:i]
                break

    def plot(self):
        fig = plt.figure(figsize=(14, 10))
        ax = fig.add_subplot(111, projection='3d')

        # Plot Central Body
        x, y, z = create_sphere(self.central["radius"], [0,0,0])
        ax.plot_surface(x, y, z, cmap=self.central['cmap'])

        # Plot Trajectories
        for i in range(self.n_body):
            #ax.plot(self.u[0, i], self.u[0, i + self.n_body], self.u[0, i + 2 * self.n_body], marker="X", color='k')
            ax.plot(self.u[:, i], self.u[:, i + self.n_body], self.u[:, i + 2 * self.n_body],
                    label=f"Trajectory {self.bodies[i]['name']}")

            if "radius" in self.bodies[i]:
                R = self.bodies[i]['radius']
                coord = [self.u[-1, i], self.u[-1, i + self.n_body], self.u[-1, i + 2 * self.n_body]]
                x, y, z = create_sphere(R, coord)
                ax.plot_surface(x,y,z, cmap=self.bodies[i]['cmap'])

        lim = np.max(np.abs(self.u[:, :3]))
        ax.set_xlim([-lim, lim])
        ax.set_ylim([-lim, lim])
        ax.set_zlim([-lim, lim])

        ax.set_xlabel('x-axis [m]')
        ax.set_ylabel('y-axis [m]')
        ax.set_zlabel('z-axis [m]')

        ax.legend()
        plt.show()

    def plot_orbital_params(self, sample_interval, format='hour'):
        sample_step = int(sample_interval / self.dt)
        size = (self.N - 1) // (sample_step + 1) + 1
        koe = np.zeros((size, 7))
        for i,j in enumerate(range(0, self.N, sample_step)):
            cur = np.array(CSV_TO_KOE(self.central['mass'], self.u[j]))
            koe[i,:6] = cur
            koe[i, 6] = i * sample_interval

        koe[:, 0] = (koe[:, 0] - self.central['radius']) / 1e3 # convert semi-major axis to altitude in km

        if format == 'hour':
            koe[:, 6] /= 3600
        elif format == 'day':
            koe[:, 6] /= 3600 / 24

        fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2,3, sharex='col')
        days = self.T_tot / 24 / 60 / 60
        fig.suptitle(f'Keplerian Orbital Elements of {self.bodies[0]["name"]} over a period of {days} days')

        #ax1.ticklabel_format(useOffset='False')
        ax1.plot(koe[:, 6], koe[:, 0])
        ax1.set_title('Semi-major axis [km]')
        ax1.grid(True)

        #ax2.ticklabel_format(useOffset='False')
        ax2.plot(koe[:, 6], koe[:, 1])
        ax2.set_title('Eccentricity [-]')
        ax2.grid(True)

        #ax3.ticklabel_format(useOffset='False')
        ax3.plot(koe[:, 6], koe[:, 2])
        ax3.set_title('Inclination [$^\circ$]')
        ax3.grid(True)

        #ax4.ticklabel_format(useOffset='False')
        ax4.plot(koe[:, 6], koe[:, 3])
        ax4.set_title('Argument of periapsis [$^\circ$]')
        ax4.set_xlabel('Time [h]')
        ax4.grid(True)

        #ax5.ticklabel_format(useOffset='False')
        ax5.plot(koe[:, 6], koe[:, 4])
        ax5.set_title('Longitude of ascending node [$^\circ$]')
        ax5.set_xlabel('Time [h]')
        ax5.grid(True)
        
        #ax6.ticklabel_format(useOffset='False')
        ax6.plot(koe[:, 6], koe[:, 5])
        ax6.set_title('Mean anomaly [$^\circ$]')
        ax6.set_xlabel('Time [h]')
        ax6.grid(True)


        plt.show()
