eps0 = 8.854187817e-12

a = 0.016
b = 0.013
h = a-b
pi = 3.14159265359
min_gap = 1e-6
d = 40e-6

stiff = 10
elnum = 8

V_amp = 100
V_freq = 62000
res = 400
ind = 20e-3

V_amp = 1000
V_freq = 500e3
res = 3636
ind = 0.08


phd1 = 0
phd2 = -1/V_freq/6
phd3 = -1/V_freq/3

phd1 = 0
phd2 = pi/3
phd3 = -pi/3

mass = 0.00001
g = 9.8


dt = 1/V_freq/16
dt_hold = dt
time_hold = 1/V_freq*32
time = time_hold  + 1/V_freq*4048

S = 2*3.1415*a*h
C0 = eps0*S