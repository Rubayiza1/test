import matplotlib.pyplot as plt
import math
import numpy as np
import black
import scipy

# Simulating two body motion under gravitational forces

# init variables
G = 6.67*10**(-11)
m1 = 5.972*10**(24)
m2 = 7.35*10**(22)
x = 385000*10**3  # (m)
x_pos = []
y = 0  # (m)
y_pos = []
Vx = 0  # (m/s)
x_vel = []
Vy = 1080  # (m/s)
y_vel = []
t = 0  # (s)
time = []
theta = math.atan2(y,x)
O = []
Ax = []
Ay = []
dt = (86400/1000)
print("theta start" = theta)

# time-evolution 'explicit' approximation
while t < 2919999:
    x_pos.append(x)
    y_pos.append(y)
    x_vel.append(Vx)
    y_vel.append(Vy)
    time.append(t)
    Ax.append(ax)
    Ay.append(ay)
    r = math.sqrt((x**2)+(y**2))
    F = -G*m1*m2/r**2
    a = F/m2
    ax = a*math.cos(theta)
    ay = a*math.sin(theta)
    x = x + Vx*dt
    y = y + Vy*dt
    Vx = Vx + ax*dt
    Vy = Vy + ay*dt
    theta = math.atan2(y,x)
    O.append(theta)
    t = t + dt


# time-evolution 'implicit' approximation !!doesn't work!!
while t < 3293200:
    x_pos.append(x)
    y_pos.append(y)
    x_vel.append(Vx)
    y_vel.append(Vy)
    time.append(t)
    r = math.sqrt(x**2+y**2)
    theta = math.atan2(y,x)
    O.append(theta)
    F = G*m1*m2/r**2
    a = F/m2
    ax = a*math.cos(theta)
    ay = a*math.sin(theta)
    Vx = Vx + ax*dt
    Vy = Vy + ay*dt
    x = x + Vx*dt
    y = y + Vy*dt
    t = t + dt

# Plotting Variables

plt.figure(1)
plt.plot(x_pos,y_pos)
plt.title("Orbit")
plt.xlabel("X-position")
plt.ylabel("Y-position")
plt.show()


plt.figure(2)
plt.plot(time, PE, time, KE, time, TE)
plt.title("Energy vs Time")
plt.xlabel("Time")
plt.ylabel("Energy")
plt.legend(["Potental", "Kinetic", "Total"], loc="upper right")
plt.show()

# Numerical Integration

# Finding the total energy of the system

#Calculating radial distance 
r_pos = np.sqrt(np.square(x_pos) + np.square(y_pos))

# Formulae
PE = -G*m1*m2/r_pos
KE = G*m1*m2/(2*r_pos)
TE = PE + KE

# Forward Reimann Sum
i = 0
U = 0
T = 0
H = 0
while i < len(r_pos):
    u = PE[i]*dt
    U = U + u
    t = KE[i]*dt
    T = T + t
    h = TE[i]*dt
    H = H + h
    i += 1

print("<<Reimann Sum>>")
print("Total Kinetic E = ",T, "Total Potential E = ",U, "Total Energy =", H)
print("check T + V - U =", T+U-H)


# Trapezoidal Sum
j = 0
U = 0
T = 0
H = 0
while j+1 < len(r_pos):
    u = PE[j]*dt + 0.5*dt*(PE[j]-PE[j+1])
    U = U + u
    t = KE[j]*dt + 0.5*dt*(KE[j]-KE[j+1])
    T = T + t
    h = TE[j]*dt + 0.5*dt*(TE[j]-TE[j+1])
    H = H + h
    j += 1

print("<<TRAPEZOIDAL METHOD>>")
print("Total Kinetic E = ",T, "Total Potential E = ",U, "Total Energy =", H)
print("check T + V - U =", T+U-H)


    
