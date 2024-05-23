#NBVAL_IGNORE_OUTPUT
from devito import Grid, TimeFunction, Eq, solve, configuration, Operator, Constant, Eq, solve, mmax
from sympy.abc import a
from sympy import simplify, sin
import numpy as np
import time
from math import pi, exp

configuration['safe-math'] = 1

# Some variable declarations
Nx = 4000
Ny = 4000
T = 1.0
nt = 1000
dt = 0.001
mu = 0.1

dx = 1.0
dy = 1.0
Lx = dx*(Nx-1)
Ly = dy*(Ny-1)

# Initialize `u` for space order 2
grid = Grid(shape=(Nx, Ny), extent=(Lx, Ly), dtype = np.float64)
u = TimeFunction(name='u', grid=grid, space_order=2)

#initial condition
x, y = grid.dimensions
init_eq = Eq(u, 5*sin(pi*(x)/Lx)*sin(pi*(y)/Ly), subdomain=grid.interior)
init_eq
init_op = Operator(init_eq, opt='noop')
init_op(time=1)
#print(u.data[0])


## Create an operator with second-order derivatives
a = Constant(name='a')
eq = Eq(u.dt, a * u.laplace, subdomain=grid.interior)
stencil = solve(eq, u.forward)
eq_stencil = Eq(u.forward, stencil)

# Create boundary condition expressions
t = grid.stepping_dim

bc = [Eq(u[t+1, 0, y], 0.)]  # left
bc += [Eq(u[t+1, Nx-1, y], 0.)]  # right
bc += [Eq(u[t+1, x, Ny-1], 0.)]  # top
bc += [Eq(u[t+1, x, 0], 0.)]  # bottom

#main finite difference operator
#op = Operator([eq_stencil]+bc, platform='nvidiaX', opt=('advanced', {'gpu-fit': u}))
op = Operator([eq_stencil]+bc, opt=('advanced', {'openmp': True}))

start_time = time.time() 

op(time=nt, dt=dt, a=mu)

end_time = time.time()

elapsed_time = end_time - start_time

print(f"Runtime: {elapsed_time} seconds")

#Solution checking
u_analytical = TimeFunction(name='u_analytical', grid=grid, space_order=2, dtype=np.float64)
error_field = TimeFunction(name='error_field', grid=grid, space_order=2, dtype=np.float64)
u_analytical = 5*exp(-mu*pi*pi*(1/Lx/Lx+1/Ly/Ly)*nt*dt)*sin(pi/Lx*x)*sin(pi/Ly*y)

error_eq = Eq(error_field, abs((u_analytical-u)/u_analytical), subdomain=grid.interior)

check_op = Operator(error_eq, opt='noop')
check_op(time=1)
print(mmax(error_field))

