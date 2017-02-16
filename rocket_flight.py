"""
Sathyan Padmanabhan, MS 
#1 Rocket Flight
Solving equation of motion to determine the rocket flight
Used second order Runge-Kutta method for solving ODE
The assignment is part of Practical Numerical Methods with Python (MAE6286), George Washington University
Program prints the flight path of the rocket in terminal, we can save the results by printing > to file
https://github.com/Sathyan64 

"""

import numpy 
from matplotlib import pyplot

mp0 = 100 #initial weight of the rocket propellant
T = 40 # time in seconds
dt = 0.001 #timestep size
N = int(T/dt)+1 #Number of nodes

#initial conditions
v0 = 0.0 
h0 = 0.0
mdotp0 = 20.0 
g = 9.8
ms = 50 # weight of the rocket shell
rho = 1.091 #average air density
A = numpy.pi*0.5**2 #maximum cross section area of the rocket
ve = 325 #exhaust velocity in m/s
cd = 0.15 #drag coefficient

u = numpy.empty((N, 2)) #initializing array size
u[0] = numpy.array([h0,v0]) #initial conditions

def f(u,mp,mdotp):
    v = u[1] 
    return numpy.array([v, (-g+ ((mdotp*ve)/(ms+mp)) - ((0.50*rho*v*numpy.abs(v)*A*cd)/(ms+mp)))])

def rkutta(u,dt,f,mp,mdotp):
    """Returns the solution at the next time-step using Runge-Kutta method.
    
    Parameters
    ----------
    u : array of float
        solution at the previous time-step.
    f : function
        function to compute the right hand-side of the system of equations.
    dt : float
        time-increment.
    u_star : Takes the midpoint between the time step size and the approximate solution 
    at the time step

    Returns
    -------
    u_n_plus_1 : array of float
        approximate solution at the next time step.
    """
    u_star = u + 0.5*dt*f(u,mp,mdotp) #midpoint method
    return u + dt * f(u_star,mp,mdotp)

for n in range (0,N-1):
    if n < (5+dt)/dt: # variable burning rate function
        mdotp = 20
        mp = mp0 - 20*n*dt 
    else:
        mp = 0
        mdotp = 0
    u[n+1] = rkutta(u[n],dt,f,mp,mdotp) #calls RK method to solve for the next step


H = u[:,0] # storing the height variable
V = u[:,1] # storing the velocity variable


for n in range (0,N-1):
    print ("H = %.4f, V = %.4f, Time = %.4f" %(H[n],V[n],n*dt))  #prints Height, velocity for the corresponding time


 



         
