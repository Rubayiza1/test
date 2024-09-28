import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.integrate import quad
import black

# Simulating two body motion under gravitational forces

# init variables
G = 6.67*10**(-11)  # Universal Gravitational Constant
m1 = 5.972*10**(24)  # mass of earth
m2 = 7.35*10**(22)  # mass of moon
x = 385000*10**3  # initial x-position of moon
x_pos = []
y = 0  # input initial y-position of moon
y_pos = []
Vx = 0  # input initial x-velocity of moon
x_vel = []
Vy = 1080  # input initial y-velocity of moon
y_vel = []
t = 0  # (s)
time = []
theta = math.atan2(y,x)
O = []
Ax = []
Ay = []
dt = (86400/1000)  #  input delta t (tune precision) - lower dt will take longer
print("theta start = ",theta)

# time-evolution 'explicit' approximation
while t < 2919999:
    x_pos.append(x)
    y_pos.append(y)
    x_vel.append(Vx)
    y_vel.append(Vy)
    time.append(t)
    r = math.sqrt((x**2)+(y**2))
    F = -G*m1*m2/r**2
    a = F/m2
    ax = a*math.cos(theta)
    ay = a*math.sin(theta)
    Ax.append(ax)
    Ay.append(ay)
    x = x + Vx*dt
    y = y + Vy*dt
    Vx = Vx + ax*dt
    Vy = Vy + ay*dt
    theta = math.atan2(y,x)
    O.append(theta)
    t = t + dt

# plots
plt.figure(1)
plt.plot(x_pos,y_pos)
plt.title("Orbit (Euler explicit approximation)")
plt.xlabel("X-position")
plt.ylabel("Y-position")
plt.show()

# time-evolution rk4 approximation 

# rk4 approximation of an elliptical orbit given initial conditions and elliptical dimensions a and b.

# finding major and minor axes
a = max(x_pos)-min(x_pos)  # semi-major axis
b = max(y_pos)-min(y_pos)  # semi-minor axis
dx = -385000000/100000  # input delta-x (precision) - lower will take longer

def dydx(x, y, a, b):
    return (x/y)*(1-(b**2/a**2))


def rk4(x0, y0, a, b, dx):
    x_pos = []
    y_pos = []
    a = a
    b = b
    x = x0
    y = y0
    theta = 0
    while x > 0:
        k1 = dx*dydx(x, y, a, b)
        k2 = dx*dydx(x + 0.5*dx, y + 0.5*k1, a, b)
        k3 = dx*dydx(x + 0.5*dx, y + 0.5*k2, a, b)
        k4 = dx*dydx(x + dx, y + k3, a, b)
        y = y + (1/6)*(k1+2*k2+2*k3+k4)
        y_pos.append(y)
        x = x + dx
        x_pos.append(x)
    return 

rk4(385000000, 1, a, b, dx)

# plots
plt.figure(2)
plt.plot(x_pos,y_pos)
plt.title("Orbit (RK4 approximation)")
plt.xlabel("X-position")
plt.ylabel("Y-position")
plt.show()

print(a,b)

# Numerical Integration

# Finding the total energy of the system

#Calculating radial distance 
r_pos = np.sqrt(np.square(x_pos) + np.square(y_pos))

# Formulae
PE = -G*m1*m2/r_pos
KE = G*m1*m2/(2*r_pos)
TE = PE + KE

# Forward Reimann Sum
R_sum = []
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
    R_sum.append(H)
    i += 1

print("<<Reimann Sum>>")
print("Kinetic E = ",T, "Potential E = ",U, "Total E=", H)
print("check T + V - U =", T+U-H)


# Trapezoidal Sum
T_sum = [0]
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
    T_sum.append(H)
    j += 1

print("<<TRAPEZOIDAL METHOD>>")
print("Kinetic E = ",T, "Potential E = ",U, "Total E =", H)
print("check T + V - U =", T+U-H)

# Simpson sum
S_sum = [0]
k = 0
U = 0
T = 0
H = 0
while k+1 < len(r_pos):
    u = (dt/6)*(PE[k] + 2*(PE[k]+PE[k+1]) + PE[k+1])
    U = U + u
    t = (dt/6)*(KE[k] + 2*(KE[k]+KE[k+1]) + KE[k+1])
    T = T + t
    h = (dt/6)*(TE[k] + 2*(TE[k]+TE[k+1]) + TE[k+1])
    H = H + h
    S_sum.append(H)
    k += 1
    
print("<<SIMPSON's METHOD>>")
print("Kinetic E = ",T, " Potential E = ",U, "Total E =", H)
print("check T + V - U =", T+U-H)

plt.figure(3)
plt.plot(time, R_sum, time, T_sum, time, S_sum)
plt.title("Riemann, Trapezoidal and Simpson Sums for Total Energy")
plt.xlabel("time")
plt.ylabel("Energy")
plt.legend(["Riemann","Trapezoidal","Simpson"],loc="upper right")
plt.show()

print("Curves overlapping may not show")

# Scipy integration of total energy

def integrand(r):
    return -G*m1*m2/(r)

TE,err = quad(integrand, min(r_pos), max(r_pos))
print("Total Energy w SciPy = ",TE)
print("There does seem to be a significant discrepancy between my code and that of SciPy")
print("Not entirely sure if this a nature of the algorithm or a fundamental flaw in my code :(")
