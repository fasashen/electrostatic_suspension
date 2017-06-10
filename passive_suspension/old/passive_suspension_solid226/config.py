l = 0.04
h = 100e-6
w = 0.04
S = w*l
elsize = 0.02
eps_0 = 8.854e-12

V_amp = 100
V_freq = 300000
V_delay = -1/V_freq/4
V_ic = 2

res = 100
ind = 100e-3

tolerance = 1e-3

dt = 1/V_freq/8
time = 1/V_freq*5000