from numpy import pi, sin, sinh, cos, arange, subtract, around, exp, insert,arctan
import matplotlib.pyplot as plt
from config import *

def calc_theory_force(gap_list):
    force = []
    for gap in gap_list:
        C = EPS_0*S/gap
        peak_coef = V_freq*ind-1/V_freq/C
        print(V_freq*ind, 1/V_freq/C) 
        q0 = V_amp/(V_freq*(res**2+(peak_coef)**2)**0.5)
        force.append( q0**2/(2*EPS_0*S) )
    return force
def problem_properties():
    gap = 10e-6
    C = EPS_0*S/gap
    phi0 = arctan(res/(V_freq*ind-1/V_freq/C))
    q0 = V_amp/(V_freq*(res**2+(V_freq*ind-1/V_freq/C)**2)**0.5)
    q0_2 = V_amp*S*EPS_0*4*pi/(  (ind*S*EPS_0*4*pi*V_freq**2-4*pi*gap)**2 + (res*S*EPS_0*4*pi*V_freq)**2  )**0.5
    print('Problem properties with gap:',gap)
    print('Capacitance:',C,'| R:',res,'| L:',ind,'| w:',V_freq, '| wL-1/wc',V_freq*ind -1/V_freq/C, '| phi0 =',phi0,'| sin(phi0)=',sin(phi0), '| xi =',1/2*res*(C/ind)**0.5)
    print('q0 =',q0, '| U_cmax =',q0/C)


gap_list = list(arange(1e-6,20e-6,0.1e-6))
force_theory = calc_theory_force(gap_list)
problem_properties()
fig = plt.figure()
fig.set_size_inches(6, 3)
p1 = fig.add_subplot(111)
p1.plot(gap_list,force_theory,marker='o', label="",linewidth=1,linestyle='--', color='b',markersize=4, mew=0)
fname = 'plottest'
fig.savefig(fname)
