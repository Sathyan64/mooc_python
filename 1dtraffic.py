"""
Sathyan Padmanabhan, MS 
#2 1dtraffic.py
Solving 1Dtraffic PDE equations using finite difference method. Time is discretized using forward and the space 
using backward difference methods.
The assignment is part of Practical Numerical Methods with Python (MAE6286), George Washington University
Program also prints the velocity of the cars in terminal, we can save the results by printing > to file
https://github.com/Sathyan64 
"""

import numpy
from matplotlib import pyplot



L = 11 #km
Vmax = 136/60 #km/min
rhomax = 250 #cars/km
nx = 51
dx = L/(nx-1)
dt = 0.001*60 #min
nt = int(3.0/dt) 

#initial condition declaration
X = numpy.linspace(0,L,nx)
rho = numpy.ones(nx)*20
rho[10:20] = 50

V = numpy.ones(nx)

"""
∂ρ/∂t+∂F/∂x = 0
F is the flux
∂F/∂x can be expanded as ∂F/∂ρ*∂ρ/∂x (chain rule)
∂F/∂ρ can be calculaed from the expression 
F(ρ)=Vmax*ρ(1−ρ/ρmax).
Sympy function is used separately to calculate ∂F/∂ρ.
∂ρ/∂t+(∂F/∂ρ*∂ρ/∂x) = 0 is discretized below 
"""
for n in range(1,nt):  
    rhon = rho.copy()
    rho[1:] = rhon[1:]-(Vmax- 2.*rhon[1:]*Vmax/rhomax)*dt/dx*(rhon[1:]-rhon[0:-1]) 
    rho[0] = 20.0

print ("X \t V(m/s)")

"""
Velocity of the cars is expressed in terms of ρ as
V = Vmax*(1−ρ/ρmax). (Km/hr)
To convert the V in (m/s) multiply the above equation 
by (1000/60)
"""

for i in range(0,nx):
    V[i] = Vmax*(1 - (rho[i]/rhomax))*(1000/60)
    print("%.4f \t %.6f" %(i*dx, V[i]))

ave = numpy.average(V)

print("Average at time : %.4f" %(ave)) 

pyplot.plot(X, V, color='#003366', ls='--', lw=3)
pyplot.show()

