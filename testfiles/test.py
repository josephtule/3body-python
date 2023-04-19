import matplotlib.pyplot as plt
import numpy as np

from gekko import GEKKO

# create GEKKO model
m = GEKKO(remote=False)
m.open_folder()
nt = 101
m.time = np.linspace(0,500,nt)

# Variables

# Initial Position and Velocity - Chaser
r1_c = -1182.339411   # [km] 
r2_c = 6816.939420    # [km] 
r3_c = 904.891745     # [km] 
v1_c = 0.1175776      # [km/s] 
v2_c = -0.963776      # [km/s] 
v3_c = 7.494102       # [km/s] 

# Initial Position and Velocity - Target
r1_t = -1182.959348    # [km] 
r2_t = 6817.396210     # [km] 
r3_t = 904.495486      # [km] 
v1_t = 0.1175776       # [km/s] 
v2_t = -0.963776       # [km/s] 
v3_t = 7.494102        # [km/s] 

# initial values
x =  m.Var(value = r1_c)
y =  m.Var(value = r2_c)
z =  m.Var(value = r3_c)
vx =  m.Var(value = v1_c)
vy =  m.Var(value = v2_c)
vz =  m.Var(value = v3_c)


# initial values
x2 =  m.Var(value = r1_t)
y2 =  m.Var(value = r2_t)
z2 =  m.Var(value = r3_t)
vx2 =  m.Var(value = v1_t)
vy2 =  m.Var(value = v2_t)
vz2 =  m.Var(value = v3_t)

# Cost Initial - Integrated thrust (fuel usage)
J =  m.Var(value = 0)

# Manipulated Variable
U = m.MV(value=0, lb=-0.01, ub=0.01)
U.STATUS = 1

# Mask time
p = np.zeros(nt) # mask final time point
p[-1] = 500
final = m.Param(value=p)

# Parameters
RT =  m.Const(6378.139);            # [km] 
mu_t = m.Const(3.9860064*10**(5))   # [km^3/s^2] 

# Equations

# Define intermediate quantities - Chaser
Rc = m.Intermediate((x**2 + y**2 + z**2)**0.5)
v = m.Intermediate((vx**2 + vy**2 + vz**2)**0.5)

ax = m.Intermediate(x * -mu_t / Rc**3 + U *vx / v )
ay = m.Intermediate(y * -mu_t / Rc**3 + U *vy / v )
az = m.Intermediate(z * -mu_t / Rc**3 + U *vz / v )

# Define intermediate quantities - Target
Rt = m.Intermediate((x2**2 + y2**2 + z2**2)**0.5)
v2 = m.Intermediate((vx2**2 + vy2**2 + vz2**2)**0.5)

ax2 = m.Intermediate(x2 * -mu_t / Rt**3 )
ay2 = m.Intermediate(y2 * -mu_t / Rt**3 )
az2 = m.Intermediate(z2 * -mu_t / Rt**3 )


# Governing equations
m.Equations((vx.dt() == ax, vy.dt() == ay, vz.dt() == az))
m.Equations((x.dt() == vx, y.dt() == vy, z.dt() == vz))

m.Equations((vx2.dt() == ax2, vy2.dt() == ay2, vz2.dt() == az2))
m.Equations((x2.dt() ==  vx2, y2.dt() ==  vy2, z2.dt() == vz2))

# Equation relating thrust to fuel usage
#m.Equation(J.dt() == m.abs2(U))
m.Equation(J.dt()**2 == U**2)


# Path Constraints

# specify endpoint conditions
# soft constraints
# m.Minimize(final*(v-v2)**2)
# m.Minimize(final*(x-x2)**2)
# m.Minimize(final*(y-y2)**2)
# m.Minimize(final*(z-z2)**2)
# m.Minimize(final*(vx-vx2)**2)
# m.Minimize(final*(vy-vy2)**2)
# m.Minimize(final*(vz-vz2)**2)    

# Constraints 
# m.abs2(Rc - RT) > 200
# m.abs2(Rc - RT) < 1000
# m.abs2(Rc*final - Rt*final) == 0
RcRT = m.Var(lb=200,ub=1000)
m.Minimize(RcRT**2 - (Rc-RT)**2)

# Objective, minimizing fuel usage
m.Minimize(J * final)

# Set solver mode - Optimal Control
m.options.IMODE = 6

# Increase maximum number of allowed iterations
m.options.MAX_ITER = 2000

# Set number of nodes per time segment
m.options.NODES = 3

# Run solver and display progress
m.solve(disp=True)