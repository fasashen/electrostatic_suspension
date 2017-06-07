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

def start_simulation(name):
    os.system(name)

def read_results():
    time, ux, fx_pos, fx_neg = [], [], [], []
    with open('./output/master_UX.txt','r') as f:
        for line in f.readlines()[14:]:
            time.append(float(line[:15]))
            ux.append(float(line[16:]))
    with open('./output/Fx_pos.txt','r') as f:
        for line in f.readlines()[14:]:
            fx_pos.append(float(line[16:]))
    with open('./output/Fx_neg.txt','r') as f:
        for line in f.readlines()[14:]:
            fx_neg.append(float(line[16:]))
    return time, ux, fx_pos, fx_neg
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
# V_amp = V_amp/2
# plot_dynamic()

start_simulation('sphere_gyroscope_SEE3.bat')

time, ux, fx_pos, fx_neg = read_results()

fx_sum = []
for i,k in zip(fx_pos, fx_neg):
    fx_sum.append(i+k)

print(time)

fig = plt.figure()
fig.set_size_inches(12, 14)

p1 = fig.add_subplot(311)
p1.plot(time,ux,label='$ux$',linewidth=1, color='r')
p1.set_xlabel("$time, s$", fontsize=15)
p1.set_ylabel("$ux, m$", fontsize=15)
p1.legend(loc='upper right',fontsize=15,numpoints=1)
p1.grid()

p2 = fig.add_subplot(312)
p2.plot(time,fx_pos,label='$pos$',linewidth=1, color='r')
p2.plot(time,fx_neg,label='$neg$',linewidth=1, color='b')

p2.set_xlabel("$time, s$", fontsize=15)
p2.set_ylabel("$force, N$", fontsize=15)
p2.legend(loc='upper right',fontsize=15,numpoints=1)
p2.grid()

p2 = fig.add_subplot(313)
p2.plot(time,fx_sum,label='$sum$',linewidth=2, color='k')
p2.set_xlabel("$time, s$", fontsize=15)
p2.set_ylabel("$force, N$", fontsize=15)
p2.legend(loc='upper right',fontsize=15,numpoints=1)
p2.grid()

fig.savefig('plot',dpi=300)