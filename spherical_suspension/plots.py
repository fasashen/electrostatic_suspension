from numpy import pi, sin, sinh, cos, arange, subtract, around, exp, insert, arctan
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from config import *
import os
  
def read_from_file(fname,header=False):
    '''
    Reads the results of capacitance calculation
    '''
    cap_list, gap_list = [], []
    with open(fname,'r') as f:
        if header == True: f = f.readlines()[14:]
        for line in f:
            gap_list.append(float(line.split()[0]))
            cap_list.append(float(line.split()[1]))
    return cap_list, gap_list

def plot_dynamic():
    # uy_list, time_list = read_from_file('./output/dyn_UY_vs_Time.txt',True)
    force_list, time_list = read_from_file('./output/master_side201_FORCE.txt',True)
    volt_list, time_list = read_from_file('./output/master_side201_VOLT.txt',True)

    C = eps0*S/d*2
    e_ode = odeint(charge, [0.0, 0.0], time_list)[:, 0]
    U_ode = [x/C for x in e_ode]
    F_ode = [x**2/(2*eps0*S)/609 for x in e_ode]

    fig = plt.figure()
    fig.set_size_inches(12, 8)
    # p1 = fig.add_subplot(311)
    # p1.plot(time_list ,[-init_gap+x for x in uy_list],marker='o', label="FEM trans126",linewidth=2, color='r',markersize=0, mew=0)
    # p1.plot([0, time_list[-1]],[0,0],linewidth=2, color='k')
    
    # p1.set_xlabel("time, sec", fontsize=15)
    # p1.set_ylabel("displacement, m", fontsize=15)


    p2 = fig.add_subplot(211)
    p2.plot(time_list ,force_list,marker='o', label="CAE TRANS126",linewidth=2, color='r',markersize=0, mew=0)
    p2.plot(time_list ,F_ode,marker='o', label="Passive susp. analogue",linewidth=0.5, color='b',markersize=0, mew=0,linestyle='-')
    p2.legend(loc='best',fontsize=15,numpoints=1)
    p2.grid()
    p2.set_xlabel("time, sec", fontsize=15)
    p2.set_ylabel("forces sum, N", fontsize=15)
    p3 = fig.add_subplot(212)
    p3.grid()
    
    p3.plot(time_list ,volt_list,marker='o', label="CAE TRANS126",linewidth=2, color='r',markersize=0, mew=0)
    p3.plot(time_list ,U_ode,marker='o', label="Passive susp. analogue",linewidth=0.5, color='b',markersize=0, mew=0,linestyle='-')
    p3.legend(loc='best',fontsize=15,numpoints=1)
    p3.set_xlabel("time, sec", fontsize=15)
    p3.set_ylabel("voltage drop on capacitor, V", fontsize=15)
    fig.savefig('dynamic_plot2',dpi=300)

def charge(y, t):
     y1, y2 = y
     dydt = [y2, 1/ind*(V_amp*sin(2*pi*V_freq*t) - y1*d/S/eps0 - res*y2)]
     return dydt

'''SOLVE CAPACITY'''

# cap_cae, gap_list = calc_cap_vs_gap(start_gap, end_gap, n)

# problem_properties()

# '''SOLVE FORCE'''
# force_cae, gap_list = calc_force_vs_gap(start_gap, end_gap, n)
# save_solution(gap_list,cap_cae,force_cae)

# plot_cap()
# plot_force_vs_gap()

# mass =  0.122395/2/g
# init_gap =  3e-6
     
# uy_list, time_list = calc_dynamic(init_gap,mass,dt,dt*2)
V_amp = V_amp/2
plot_dynamic()