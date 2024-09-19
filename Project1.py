
import matplotlib.pyplot as plt

# Initializing variables

m1 = 100
m2 = 1
rp = 10000
vp = 7000
G = 1
r = rp
v = vp
dt = 10

while t < 2000000:  # Trajectory
    F = G*m1*m2/r
    a = F/m2
    v = v + a*dt
    r = r + v*dt
    t = t + dt
    
