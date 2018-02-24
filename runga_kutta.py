# https://rosettacode.org/wiki/Runge-Kutta_method#Python

def RK4(f):
    return lambda t, y, dt: (
        lambda dy1: (
            lambda dy2: (
                lambda dy3: (
                    lambda dy4: (dy1 + 2 * dy2 + 2 * dy3 + dy4) / 6
                )(dt * f(t + dt, y + dy3))
            )(dt * f(t + dt / 2, y + dy2 / 2))
        )(dt * f(t + dt / 2, y + dy1 / 2))
    )(dt * f(t, y))


def theory(t): return (t ** 2 + 4) ** 2 / 16


from math import sqrt

dy = RK4(lambda t, y: t * sqrt(y))

t, y, dt = 0., 1., .1
while t <= 10:
    if abs(round(t) - t) < 1e-5:
        print("y(%2.1f)\t= %4.6f \t error: %4.6g" % (t, y, abs(y - theory(t))))
    t, y = t + dt, y + dy(t, y, dt)

# possible alternative solution https://github.com/MechCoder/Rk4-Polynomials/blob/master/rk4.py

import math

m1 = 0.2
m1 = m2
L1 = 6*0.0254
L2 = L1
l1 = L1/2
l2 = L2/2
g = 9.81

Tao = 0.1*sin(5*t)
b1 = 0.001
b2 = 0.001
J1 = m1*L1^(2/3)
J2 = m2*L2^(3/3)

J0_h = J1 + m1*l1^2 + m2*L1^2
J1_h = J1 + m1*l1^2
J2_h = J2 + m2*l2^2

t1 =
t2 =
t1_d =
t2_d =

# B_dd EQNS
F = (-sin(2*t2)*J2_h^(2*t1_d*t2_d)) - (
    L1*(sin(t2) - sin(t2)^3)*J2_h*l2*m2*t1_d^2) + (
    L1*sin(t2)*J2_h*12*m2*t2_d^2) - (
    b1-J2_h*t1_d) + (
    Tao*J2_h) + (
    L1*g*sin(2*t2)*l2^2*m2^(2/2)) + (
    L1*b2*cos(t2)*l2*m2*t2_d)

G = J2_h^2*sin(t2)^2 + J0_h