import numpy 
from matplotlib import pyplot
from math import sin, cos, log, ceil

dt_values = numpy.array([0.1, 0.05, 0.01, 0.005, 0.001, 0.0001])
# array that will contain solution of each grid
u_values = numpy.empty_like(dt_values, dtype=numpy.ndarray)

for i, dt in enumerate(dt_values):
    mp0 = 100
    T = 40
    N = int(T/dt)+1
    v0 = 0.0
    h0 = 0.0
    mdotp0 = 20.0
    g = 9.8
    ms = 50
    rho = 1.091
    A = numpy.pi*0.5**2
    ve = 325
    cd = 0.15

    u = numpy.empty((N, 2))
    u[0] = numpy.array([h0,v0])

    def f(u,mp,mdotp):
        v = u[1] 
        return numpy.array([v, (-g+ ((mdotp*ve)/(ms+mp)) - ((0.50*rho*v*numpy.abs(v)*A*cd)/(ms+mp)))])

    def rkutta(u,dt,f,mp,mdotp):
        u_star = u + 0.5*dt*f(u,mp,mdotp) #midpoint method
        return u + dt * f(u_star,mp,mdotp)

    for n in range (0,N-1):
        if n < (5+dt)/dt:
            mdotp = 20
            mp = mp0 - 20*n*dt
        else:
            mp = 0
            mdotp = 0
        u[n+1] = rkutta(u[n],dt,f,mp,mdotp)
        
    u_values[i] = u

#for n in range (0,N-1):
    #print ("H = %.4f, V = %.4f, Time = %.4f" %(H[n],V[n],n*dt)) 
  
def get_diffgrid(u_current, u_fine, dt):
    """Returns the difference between one grid and the fine one using L-1 norm.
    
    Parameters
    ----------
    u_current : array of float
        solution on the current grid.
    u_finest : array of float
        solution on the fine grid.
    dt : float
        time-increment on the current grid.
    
    Returns
    -------
    diffgrid : float
        difference computed in the L-1 norm.
    """
    
    N_current = len(u_current[:,0])
    N_fine = len(u_fine[:,0])
   
    grid_size_ratio = ceil(N_fine/N_current)
    
    diffgrid = dt * numpy.sum( numpy.abs(\
            u_current[:,1]- u_fine[::grid_size_ratio,1])) 
    
    return diffgrid

diffgrid = numpy.empty_like(dt_values)

for i, dt in enumerate(dt_values):
    print('dt = {}'.format(dt))

    ### call the function get_diffgrid() ###
    diffgrid[i] = get_diffgrid(u_values[i], u_values[-1], dt)

# log-log plot of the grid differences
pyplot.figure(figsize=(6,6))
pyplot.grid(True)
pyplot.xlabel('$\Delta t$', fontsize=18)
pyplot.ylabel('$L_1$-norm of the grid differences', fontsize=18)
pyplot.axis('equal')
pyplot.loglog(dt_values[:-1], diffgrid[:-1], color='k', ls='-', lw=2, marker='o')
pyplot.show();
         
