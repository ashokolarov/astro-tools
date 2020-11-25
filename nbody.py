import numpy as np
from constants import *
from plotter import * 


class Propagator():

    def __init__(self, bodies, N, dt):
        self.n_body = len(bodies)
        self.bodies = bodies
        self.masses = [elem["mass"] for elem in bodies]
        self.names = [elem["name"] for elem in bodies]
        self.N = N
        self.dt = dt
        self.Tmax = N * dt

    def rk4(self, f, u):
        """
        Runge-Kutta explicit 4-stage scheme - 1 step.
        
          t  - time
          u  - solution at t
        Return: approximation of y at t+dt.
        """
        k1 = f(u)
        k2 = f(u + self.dt*k1/2)
        k3 = f(u + self.dt*k2/2)
        k4 = f(u + self.dt*k3) 
        return u + self.dt/6 * (k1 + 2*k2 + 2*k3 + k4)

class NBodyPropagator(Propagator):

    def __init__(self, init, bodies, N, dt):
        super().__init__(bodies,N,dt)
        self.M_Multibody = 6 * self.n_body
        self.u = np.zeros((N, self.M_Multibody), dtype=np.float64)
        self.u[0, :] = init

    def f_multibody(self, u):
        x    = u[0:self.n_body]
        y    = u[self.n_body:2*self.n_body]
        z    = u[2*self.n_body:3*self.n_body]
        xdot = u[3*self.n_body:4*self.n_body]
        ydot = u[4*self.n_body:5*self.n_body]
        zdot = u[5*self.n_body:]
        f = np.zeros((self.n_body,6))
        f[:,0] = xdot
        f[:,1] = ydot
        f[:,2] = zdot
        for i in range(self.n_body):
            for j in range(self.n_body):
                if i != j:
                    r = np.sqrt( (x[j]-x[i])**2 + (y[j]-y[i])**2 + (z[j] - z[i])**2 )
                    f[i,3] += G*self.masses[j]*(x[j]-x[i])/r**3
                    f[i,4] += G*self.masses[j]*(y[j]-y[i])/r**3
                    f[i,5] += G*self.masses[j]*(z[j]-z[i])/r**3
        return f.T.flatten()

    def propagate(self):
        for i in range(self.N - 1):
            self.u[i+1,:] = self.rk4(self.f_multibody, self.u[i])
    
    def plot(self):
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Plot initial positions
        #for i in range(self.n_body):
            #ax.plot(self.u[0,:i], self.u[0,:i+2], self.u[0,:i+4], marker="x")

        # Plot trajectories
        for i in range(self.n_body):
            ax.plot(self.u[:,i], self.u[:,i+2], self.u[:,i+4], label=f"{self.names[i]}")

        limit = np.max(np.abs(self.u[:,3]))
        ax.set_xlim([-limit, limit])
        ax.set_ylim([-limit, limit])
        ax.set_zlim([-limit, limit])

        plt.legend()
        plt.show()

class TwoBodyPropagator(Propagator):

    def __init__(self, init, bodies, N, dt):
        super().__init__(bodies, N, dt)
        self.u = np.zeros((N, 6), dtype=np.float64)
        self.u[0, :] = init

    def f_twobody(self, u):
        x = u[0]
        y = u[1]
        z = u[2]
        xdot = u[3]
        ydot = u[4]
        zdot = u[5]
        f = np.zeros(6)
        f[0] = xdot
        f[1] = ydot
        f[2] = zdot
        r = np.sqrt( x**2 + y**2 + z**2 )
        f[3] = -G*self.masses[0]*x/r**3
        f[4] = -G*self.masses[0]*y/r**3
        f[5] = -G*self.masses[0]*z/r**3
        return f

    def propagate(self):
        for i in range(self.N - 1):
            self.u[i+1,:] = self.rk4(self.f_twobody, self.u[i])
        
    def plot(self):
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        x,y,z = create_sphere(self.bodies[0]["radius"])
        ax.plot_surface(x, y, z, color='b')

        ax.plot(self.u[0,0], self.u[0,1], self.u[0,2], marker='o', label="Inital Position")
        ax.plot(self.u[:,0], self.u[:,1], self.u[:,2], label=f"Trajectory")

        lim = np.max(np.abs(self.u[:,:3]))
        ax.set_xlim([-lim, lim])
        ax.set_ylim([-lim, lim])
        ax.set_zlim([-lim, lim])

        plt.legend()
        plt.show()

















        
