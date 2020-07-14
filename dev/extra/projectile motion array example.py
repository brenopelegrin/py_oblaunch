# Import various functions meant for numerical science
import numpy as np
from math import cos,sin,pi
# Import functionality for plotting
import matplotlib.pyplot as plt

N = 10000 # Number of time steps

# Create a uniformly spaced time-array
t=np.zeros([N+1,1])
vi = 5
angulo = pi/6

# Calculate the size of a time step
dt = 0.001

# Create empty acceleration, velocity and position arrays
a = np.zeros((N+1,2))
v = np.zeros((N+1,2))
r = np.zeros((N+1,2))

# Set initial conditions
v[0] = (vi*cos(angulo),vi*sin(angulo)) # inital velocity, m/s
r[0] = (0,0)  # initial position, m

m = 0.05 # mass, kg
g = 9.81 # acceleration of gravity, m/s^2
rho = 1.3 # air density, kg/m^3
C_D = 0.45 # drag coefficient
d = 0.02 # diameter of cannonball, m
A = pi*d**2 # cross-sectional area, m^2

def F(r, v, t):
    return (0, -m*g) - 0.5*rho*C_D*A*abs(v)*v

# Solving equations of motion iteratively
i=0
while r[i,1]+d/2 >= 0 and i<N:
    t[i+1] = t[i] + dt
    a[i] = F(r[i], v[i], t[i+1])/m
    v[i+1] = v[i] + a[i]*dt
    r[i+1] = r[i] + v[i]*dt
    i+=1

# Extract x and y coordinates
v = v[:i]
x = r[:i,0]
y = r[:i,1]
t = t[:i]

# Plot figure
plt.plot(t,y)

# Prettify the plot
plt.xlabel('Horizontal distance, [m]')
plt.ylabel('Vertical distance, [m]')
plt.title('Trajectory of a fired cannonball')
plt.grid()

# Makes the plot appear on the screen
plt.show()