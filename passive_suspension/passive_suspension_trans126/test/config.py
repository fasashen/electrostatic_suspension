EPS_0 = 8.854187817e-12
EPS_R = 1

a = 0.005
w = 1
S = w*a
b = 25e-6

start_gap = 1e-6
end_gap = 10e-6

n = 20

numEperLength = 1000
numEperGap = 4

V_amp = 10
V_freq = 40000

res = 300
ind = 20e-4

dt_hold = 1/V_freq/32
time_hold = 1/V_freq*8
time_hold = 0

dt = 1/V_freq/128
time = time_hold + 1/V_freq*64

g = 9.8

mass =  0.122395/2/g
init_gap =  6.039275029837934e-06/2

C0 = 4.425744022121485784379828e-14
C0 = EPS_0*S