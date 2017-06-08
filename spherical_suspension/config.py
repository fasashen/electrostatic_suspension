eps0 = 8.854187817e-12

a = 0.016
b = 0.013
h = a-b
pi = 3.14159265359
min_gap = 1e-6
d = 40e-6

stiff = 10
elnum = 8

V_amp = 30000
V_freq = 500e3
res = 3636
ind = 0.08

phd1 = 0
phd2 = pi/3
phd3 = -pi/3

mass = 0.000001
g = 9.8

dt = 1/V_freq/16
dt_hold = dt
time_hold = 1/V_freq*32
time = time_hold  + 1/V_freq*5120

S = 2*3.1415*a*h
C0 = eps0*S