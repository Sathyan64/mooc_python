"""
Sathyan Padmanabhan, MS 
#3 shocktube.py
Solving SOD's shock tube using finite difference method. Richtmyer method is used to solve the discretized 
1D euler equations
The assignment is part of Practical Numerical Methods with Python (MAE6286), George Washington University
Program also prints the velocity, pressure and density at the user given time   > to save in file
https://github.com/Sathyan64 
"""

import numpy
from matplotlib import pyplot
#from matplotlib import animation 


#Computing flux vector F

def computeF(u,gamma):
    u1 = u[0,:]
    u2 = u[1,:]
    u3 = u[2,:]
    return numpy.array([u2,((u2**2)/u1)+(gamma-1)*(u3-0.5*((u2**2)/u1)),(u2/u1)*(u3+(gamma-1)*(u3-0.5*((u2**2)/u1)))])

# similar to maccormack two-step method, Richtmyer is also a two step method. Strong oscillations 
# persist in this method. Give it a try :)  

def richtmyer(U,nt,dt,dx,gamma):
    Ustar = U.copy()
    Un = numpy.empty_like(U)
    for n in range(1,nt):
        F = computeF(U,gamma)
        Ustar[:,:-1] = 0.5*(U[:,1:]+U[:,:-1]) - 0.5*(dt/dx)*(F[:,1:] - F[:,:-1])
        Fstar = computeF(Ustar,gamma)
        Un[:,1:] = U[:,1:] - dt/dx*(Fstar[:,1:] - Fstar[:,:-1])
        Un[:,0] = U[:,0]
        U[:,:] = Un.copy()
    return Un  

 
    
nx = 81
dx = .25
dt = .0002   
T = 0.01
nt = int(T/dt)+1 

gamma = 1.4
x = numpy.linspace(-10,10,nx)
#print(x)
rwave = numpy.where(x >= 0.0)
lwave = numpy.where(x < 0.0)

u = numpy.ones(nx)
rho = numpy.ones(nx)
P = numpy.ones(nx)
et = numpy.ones(nx)

rho[lwave] = 1
rho[rwave] = 0.125
u[lwave] = 0
u[rwave] = 0
P[lwave] = 100*10**3
P[rwave] = 10*10**3
et = (P/((gamma-1)*rho))+(0.50*u**2)

U = numpy.array([rho,rho*u,rho*et])

un = richtmyer(U,nt,dt,dx,gamma)


density = un[0,:]
velocity = un[1,:]/un[0,:]
pressure = (gamma-1)*(un[2,:] - 0.50*(un[1,:]**2/un[0,:]))
print("x \t Velocity \t Pressure \t Density \n")
for i in range(1,nx):
    print("%.3f\t %.3f\t %.3f \t %.3f \n"%(x[i],velocity[i],pressure[i],density[i]))

