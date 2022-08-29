#! / usr / bin / python
# -*- coding: utf-8 -*-

import numpy, os, math 
from matplotlib import pyplot
import phugoid

def test00():
    # linspace
    a= numpy.linspace(4,23,43)
    print(a[5])

    # Array dimensions
    ones_array = numpy.ones( (5,17) ) 
    zeros_array = numpy.zeros( ones_array.shape )
    print(zeros_array.shape)

    # Array trigonometry
    p=7
    r = numpy.array([11.2, 4.7, 6.6])
    q = numpy.sin(p/r)**3
    print(q[1])

    # Initial values

    # Discretization

def test01():
    """
    Test about phugoid

    Parameters
    ----------
    variable : type
        desciption.

    Returns
    -------
    variable : type
        desciption.
    """
    phugoid.plot_flight_path(64.0, 16.0, 0.0)
    phugoid.plot_flight_path(64.0, 16.0, 180.0)
    phugoid.plot_flight_path(16.0, 48.0, 0.0)
    phugoid.plot_flight_path(64.0, 16.0, -90.0)

def test02():
    """
    Analysis about phugoid

    Parameters
    ----------
    variable : type
        desciption.

    Returns
    -------
    variable : type
        desciption.
    """
    # Create the time grid.
    T = 100.0  # length of the time-interval
    dt = 0.02  # time-step size
    N = int(T / dt) + 1  # number of time steps
    t = numpy.linspace(0.0, T, num=N)  # time grid

    # Set the initial conditions.
    z0 = 100.0  # altitude
    b0 = 10.0  # upward velocity resulting from gust
    zt = 100.0  # trim altitude
    g = 9.81  # acceleration due to gravity

    # Set the initial value of the numerical solution.
    u = numpy.array([z0, b0])

    # Create an array to store the elevation value at each time step.
    z = numpy.zeros(N)
    z[0] = z0

    # Temporal integration using Euler's method.
    for n in range(1, N):
        rhs = numpy.array([u[1], g * (1 - u[0] / zt)])
        u = u + dt * rhs
        z[n] = u[0]

    # Set the font family and size to use for Matplotlib figures.
    pyplot.rcParams['font.family'] = 'serif'
    pyplot.rcParams['font.size'] = 16

    z_exact = (b0 * (zt / g)**0.5 * numpy.sin((g / zt)**0.5 * t) +
           (z0 - zt) * numpy.cos((g / zt)**0.5 * t) + zt)

    # Plot the numerical solution and the exact solution.
    pyplot.figure(figsize=(9.0, 4.0))  # set the size of the figure
    pyplot.title('Elevation of the phugoid over the time')  # set the title
    pyplot.xlabel('Time [s]')  # set the x-axis label
    pyplot.ylabel('Elevation [m]')  # set the y-axis label
    pyplot.xlim(t[0], t[-1])  # set the x-axis limits
    pyplot.ylim(40.0, 160.0)  # set the y-axis limits
    pyplot.grid()  # set a background grid to improve readability
    pyplot.plot(t, z, label='Numerical', color='C0', linestyle='-', linewidth=2)
    pyplot.plot(t, z_exact, label='Analytical', color='C1', linestyle='-', linewidth=2)
    pyplot.legend();  # set the legend

    pyplot.show()    

    # Set the list of time-step sizes.
    dt_values = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0001]

    # Create an empty list that will contain the solution of each grid.
    z_values = []

    for dt in dt_values:
        N = int(T / dt) + 1  # number of time-steps
        t = numpy.linspace(0.0, T, num=N)  # time grid
        # Set the initial conditions.
        u = numpy.array([z0, b0])
        z = numpy.empty_like(t)
        z[0] = z0
        # Temporal integration using Euler's method.
        for n in range(1, N):
            rhs = numpy.array([u[1], g * (1 - u[0] / zt)])
            u = u + dt * rhs
            z[n] = u[0]  # store the elevation at time-step n+1
        z_values.append(z)  # store the elevation over the time

    # Create an empty list to store the errors on each time grid.
    error_values = []

    for z, dt in zip(z_values, dt_values):
        N = int(T / dt) + 1  # number of time-steps
        t = numpy.linspace(0.0, T, num=N)  # time grid
        # Compute the exact solution.
        z_exact = (b0 * (zt / g)**0.5 * numpy.sin((g / zt)**0.5 * t) +
                (z0 - zt) * numpy.cos((g / zt)**0.5 * t) + zt)
        # Calculate the L1-norm of the error for the present time grid.
        error_values.append(l1_error(z, z_exact, dt))

    # Plot the error versus the time-step size.
    pyplot.figure(figsize=(6.0, 6.0))
    pyplot.title('L1-norm error vs. time-step size')  # set the title
    pyplot.xlabel('$\Delta t$')  # set the x-axis label
    pyplot.ylabel('Error')  # set the y-axis label
    pyplot.grid()
    pyplot.loglog(dt_values, error_values, color='C0', linestyle='--', marker='o')  # log-log plot
    pyplot.axis('equal');  # make axes scale equally

    pyplot.show()   

def l1_error(z, z_exact, dt):
    """
    Computes and returns the error
    (between the numerical and exact solutions)
    in the L1 norm.
    
    Parameters
    ----------
    z : numpy.ndarray
        The numerical solution as an array of floats.
    z_exact : numpy.ndarray
        The analytical solution as an array of floats.
    dt : float
        The time-step size.
        
    Returns
    -------
    error: float
        L1-norm of the error with respect to the exact solution.
    """
    error = dt * numpy.sum(numpy.abs(z - z_exact))
    return error

def test03():
    # Set parameters.
    g = 9.81  # gravitational acceleration (m.s^{-2})
    vt = 30.0  # trim velocity (m.s)
    CD = 1.0 / 40  # drag coefficient
    CL = 1.0  # lift coefficient

    # Set initial conditions.
    v0 = vt  # start at the trim velocity
    theta0 = 0.0  # trajectory angle
    x0 = 0.0  # horizontal position
    y0 = 1000.0  # vertical position (altitude)

    T = 100.0  # length of the time interval
    dt = 0.1  # time-step size
    N = int(T / dt) + 1  # number of time steps

    # Create array to store the solution at each time step.
    u = numpy.empty((N, 4))
    # Set the initial conditions.
    u[0] = numpy.array([v0, theta0, x0, y0])

    # Time integration with Euler's method.
    for n in range(N - 1):
        u[n + 1] = euler_step(u[n], rhs_phugoid, dt, CL, CD, g, vt)

    # Get the glider's position over the time.
    x = u[:, 2]
    y = u[:, 3]

    # Set the font family and size to use for Matplotlib figures.
    pyplot.rcParams['font.family'] = 'serif'
    pyplot.rcParams['font.size'] = 16

    # Plot the path of the glider.
    pyplot.figure(figsize=(9.0, 4.0))
    pyplot.title('Path of the glider (flight time = {})'.format(T))
    pyplot.xlabel('x')
    pyplot.ylabel('y')
    pyplot.grid()
    pyplot.plot(x, y, color='C0', linestyle='-', linewidth=2)

    pyplot.show()

def rhs_phugoid(u, CL, CD, g, vt):
    """
    Returns the right-hand side of the phugoid system of equations.
    
    Parameters
    ----------
    u : list or numpy.ndarray
        Solution at the previous time step
        as a list or 1D array of four floats.
    CL : float
        Lift coefficient.
    CD : float
        Drag coefficient.
    g : float
        Gravitational acceleration.
    vt : float
        Trim velocity.
    
    Returns
    -------
    rhs : numpy.ndarray
        The right-hand side of the system
        as a 1D array of four floats.
    """
    v, theta, x, y = u
    rhs = numpy.array([-g * math.sin(theta) - CD / CL * g / vt**2 * v**2,
                       -g * math.cos(theta) / v + g / vt**2 * v,
                       v * math.cos(theta),
                       v * math.sin(theta)])
    return rhs

def euler_step(u, f, dt, *args):
    """
    Returns the solution at the next time step using Euler's method.
    
    Parameters
    ----------
    u : numpy.ndarray
        Solution at the previous time step
        as a 1D array of floats.
    f : function
        Function to compute the right-hand side of the system.
    dt : float
        Time-step size.
    args : tuple, optional
        Positional arguments to pass to the function f.
    
    Returns
    -------
    u_new : numpy.ndarray
        The solution at the next time step
        as a 1D array of floats.
    """
    u_new = u + dt * f(u, *args)
    return u_new

def l1_diff(u_coarse, u_fine, dt):
    """
    Returns the difference in the L1-norm between the solution on
    a coarse grid and the solution on a fine grid.
    
    Parameters
    ----------
    u_coarse : numpy.ndarray
        Solution on the coarse grid as an array of floats.
    u_fine : numpy.ndarray
        Solution on the fine grid as an array of floats.
    dt : float
        Time-step size.
    
    Returns
    -------
    diff : float
        The difference between the two solutions in the L1-norm
        scaled by the time-step size.
    """
    N_coarse = len(u_coarse)
    N_fine = len(u_fine)
    ratio = math.ceil(N_fine / N_coarse)
    diff = dt * numpy.sum(numpy.abs(u_coarse - u_fine[::ratio]))
    return diff

def test04():
    # Set parameters.
    g = 9.81  # gravitational acceleration (m.s^{-2})
    vt = 30.0  # trim velocity (m.s)
    CD = 1.0 / 40  # drag coefficient
    CL = 1.0  # lift coefficient

    # Set initial conditions.
    v0 = vt  # start at the trim velocity
    theta0 = 0.0  # trajectory angle
    x0 = 0.0  # horizontal position
    y0 = 1000.0  # vertical position (altitude)

    T = 100.0  # length of the time interval
    dt = 0.1  # time-step size
    N = int(T / dt) + 1  # number of time steps

    # Set the list of time-step sizes.
    dt_values = [0.1, 0.05, 0.01, 0.005, 0.001]

    # Create an empty list that will contain the solution of each grid.
    u_values = []

    for dt in dt_values:
        N = int(T / dt) + 1  # number of time-steps
        # Create array to store the solution at each time step.
        u = numpy.empty((N, 4))
        # Set the initial conditions.
        u[0] = numpy.array([v0, theta0, x0, y0])
        # Temporal integration using Euler's method.
        for n in range(N - 1):
            u[n + 1] = euler_step(u[n], rhs_phugoid, dt, CL, CD, g, vt)
        # Store the solution for the present time-step size
        u_values.append(u)    

    # Create an empty list to store the difference in the solution
    # between two consecutive grids.
    diff_values = []

    for i, dt in enumerate(dt_values[:-1]):
        diff = l1_diff(u_values[i][:, 2], u_values[-1][:, 2], dt)
        diff_values.append(diff)

    # Set the font family and size to use for Matplotlib figures.
    pyplot.rcParams['font.family'] = 'serif'
    pyplot.rcParams['font.size'] = 16

    # Plot the difference versus the time-step size.
    pyplot.figure(figsize=(6.0, 6.0))
    pyplot.title('L1-norm difference vs. time-step size')  # set the title
    pyplot.xlabel('$\Delta t$')  # set the x-axis label
    pyplot.ylabel('Difference')  # set the y-axis label
    pyplot.grid()
    pyplot.loglog(dt_values[:-1], diff_values,
                color='C0', linestyle='--', marker='o')  # log-log plot
    pyplot.axis('equal');  # make axes scale equally

    pyplot.show()

def test05():
    # Set parameters.
    g = 9.81  # gravitational acceleration (m.s^{-2})
    vt = 30.0  # trim velocity (m.s)
    CD = 1.0 / 40  # drag coefficient
    CL = 1.0  # lift coefficient

    # Set initial conditions.
    v0 = vt  # start at the trim velocity
    theta0 = 0.0  # trajectory angle
    x0 = 0.0  # horizontal position
    y0 = 1000.0  # vertical position (altitude)

    T = 100.0  # length of the time interval
    dt = 0.1  # time-step size
    N = int(T / dt) + 1  # number of time steps

    r = 2  # refinement ratio for the time-step size
    h = 0.001  # base grid size

    dt_values2 = [h, r * h, r**2 * h]
    u_values2 = []

    for dt in dt_values2:
        N = int(T / dt) + 1  # number of time steps
        # Create array to store the solution at each time step.
        u = numpy.empty((N, 4))
        # Set initial conditions.
        u[0] = numpy.array([v0, theta0, x0, y0])
        # Time integration using Euler's method.
        for n in range(N - 1):
            u[n + 1] = euler_step(u[n], rhs_phugoid, dt, CL, CD, g, vt)
        # Store the solution.
        u_values2.append(u)

    # Calculate f2 - f1.
    f2_f1 = l1_diff(u_values2[1][:, 2], u_values2[0][:, 2], dt_values2[1])
    # Calculate f3 - f2.
    f3_f2 = l1_diff(u_values2[2][:, 2], u_values2[1][:, 2], dt_values2[2])
    # Calculate the observed order of convergence.
    p = math.log(f3_f2 / f2_f1) / math.log(r)
    print('Observed order of convergence: p = {:.3f}'.format(p))

#________________________________________________________
# Principal
if __name__ == '__main__':
    os.system('cls')
    test05()