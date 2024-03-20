from examples.cfd import plot_field, init_hat
from examples.performance import unidiff_output, print_kernel
import numpy as np
#%matplotlib inline
import time

# Some variable declarations
nx = 4096
ny = 4096
nt = 400
nu = .5
dx = 2. / (nx - 1)
dy = 2. / (ny - 1)
sigma = .25
dt = sigma * dx * dy / nu


from devito import Grid, TimeFunction, Eq, solve
from sympy.abc import a
from sympy import nsimplify

# Initialize `u` for space order 2
grid = Grid(shape=(nx, ny), extent=(2., 2.))
u = TimeFunction(name='u', grid=grid, space_order=2)

eq = Eq(u.dt, a * u.laplace)
stencil = solve(eq, u.forward)
eq_stencil = Eq(u.forward, stencil)

eq_stencil

#NBVAL_IGNORE_OUTPUT
from devito import Operator, Constant, Eq, solve

# Reset our data field and ICs
init_hat(field=u.data[0], dx=dx, dy=dy, value=1.)

# Field initialization
u.data[0][ nx//4:nx//4+10 , ny//4:ny//4+10 ] = 2
u.data[0][ 3*nx//4:3*nx//4+10 , ny//4:ny//4+10 ] = 3
u.data[0][ nx//4:nx//4+10 , 3*ny//4:3*ny//4+10 ] = 4
u.data[0][ 3*ny//4:3*ny//4+10 , 3*ny//4:3*ny//4+10 ] = 5


# Create an operator with second-order derivatives
a = Constant(name='a')
eq = Eq(u.dt, a * u.laplace, subdomain=grid.interior)
stencil = solve(eq, u.forward)
eq_stencil = Eq(u.forward, stencil)

# Create boundary condition expressions
x, y = grid.dimensions
t = grid.stepping_dim
bc = [Eq(u[t+1, 0, y], 1.)]  # left
bc += [Eq(u[t+1, nx-1, y], 1.)]  # right
bc += [Eq(u[t+1, x, ny-1], 1.)]  # top
bc += [Eq(u[t+1, x, 0], 1.)]  # bottom

op = Operator([eq_stencil] + bc, opt=('advanced', {'openmp': True}))

#op = Operator([eq_stencil] + bc, opt=('advanced-fsg'))

#op0_b0_omp = Operator(eq, opt=('noop', {'openmp': True, 'par-dynamic-work': 100}))
#print_kernel(op)


#op = Operator([eq_stencil] + bc)
#op.apply(autotune=True)
start_time = time.time() 

op(time=2*nt, dt=dt, a=nu)
end_time = time.time()

elapsed_time = end_time - start_time

print(f"Runtime: {elapsed_time} seconds")