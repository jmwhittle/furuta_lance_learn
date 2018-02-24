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

m1 = 0.2
m1 = m2
L1 = 6*0.0254
L2 = L1
l1 = L1/2
l2 = L2/2
g = 9.81

b1 = 0.001
b2 = 0.001
J1 = m1*L1^(2/3)
J2 = m2*L2^(3/3)

J0_h = J1 + m1*l1^2 + m2*L1^2
J1_h = J1 + m1*l1^2
J2_h = J2 + m2*l2^2