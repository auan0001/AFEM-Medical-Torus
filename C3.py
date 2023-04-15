import fenics as fe
import random
import math
import matplotlib.pyplot as plt

# Mesh, function space
mesh = fe.Mesh("./2d_mesh/circle3.xml")

# FEM space
W = fe.VectorFunctionSpace(mesh, 'P', 1)

# Params
t = 0.0
T = 20
k = 0.05
r = 0.2
R = 0.5
alpha1 = 0.01
alpha2 = 0.005
rho = 2
cf = 0.024
ck = 0.055

# Lists
time = []
massloss0 = []
massloss1 = []

# IC class
class InitialConditions(fe.UserExpression):
    def eval(self, values, x):
        if abs(R-math.sqrt(x[0]**2 + x[1]**2)) <= r:
            values[0] = 1/2*rho*random.uniform(0,1)
            values[1] = rho+rho/2*random.uniform(0,1)
        else:
            values[0] = 0
            values[1] = 0

    def value_shape(self):
        return(2,)

# IC
indata = InitialConditions(degree=2)
u0 = fe.Function(W)
u0 = fe.interpolate(indata, W)

# Test, trial

u, v = fe.TrialFunction(W), fe.TestFunction(W)

# Bilinear, linear forms in a Crank-Nicholson scheme
# computing the next timestep as an average of the previous

a0 = u[0]*v[0]*fe.dx + 1/2*k*alpha1*fe.inner(fe.grad(u[0]),fe.grad(v[0]))*fe.dx + 1/2*k*cf*u[0]*v[0]*fe.dx
a1 = u[1]*v[1]*fe.dx + 1/2*k*alpha2*fe.inner(fe.grad(u[1]), fe.grad(v[1]))*fe.dx + 1/2*k*(cf+ck)*u[1]*v[1]*fe.dx

# Non-linear terms and load vector computed here
# since they are applying the previous timestep

L0 = u0[0]*v[0]*fe.dx - 1/2*k*alpha1*fe.inner(fe.grad(u0[0]),fe.grad(v[0]))*fe.dx \
    - 1/2*k*cf*u0[0]*v[0]*fe.dx - k*u0[0]*u0[1]*u0[1]*v[0]*fe.dx + k*cf*v[0]*fe.dx 

L1 = u0[1]*v[1]*fe.dx - 1/2*k*alpha2*fe.inner(fe.grad(u0[1]), fe.grad(v[1]))*fe.dx \
    - 1/2*k*(cf+ck)*u0[1]*v[1]*fe.dx + k*u0[0]*u0[1]*u0[1]*v[1]*fe.dx

a = a0 + a1
L = L0 + L1

# Output file
file = fe.File("./results/C3soln.pvd", "compressed")

# Assign IC
u = fe.Function(W)
u.assign(u0)

u_initial = fe.Function(W)
u_initial.assign(u0)

# Define an integral functional
M0 = (u_initial[0] - u[0])*fe.dx
M1 = (u_initial[1] - u[1])*fe.dx

file << (u,t)
while t<T:
    time.append(t)
    # use u0 as previous timestep
    u0.assign(u)
    A = fe.assemble(a)
    b = fe.assemble(L)
    fe.solve(A, u.vector(), b, "lu")

    # compute the functional
    mass0 = fe.assemble(M0)
    mass1 = fe.assemble(M1)
    massloss0.append(mass0)
    massloss1.append(mass1)

    file << (u,t)
    t += k

plt.title('Hormone mass loss in a 2D disc')
plt.xlabel('Time')
plt.ylabel('Mass loss')
plt.plot(time, massloss0, label='Hormone 1')
plt.plot(time, massloss1, label='Hormone 2')
plt.legend()
plt.grid(True)
plt.show()
