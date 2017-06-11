from numpy import pi, sin, sinh, cos, arange, subtract, around, exp, insert,arctan
import matplotlib
import matplotlib.pyplot as plt
from config import *
import os
from scipy.integrate import odeint


def get_force(gap_list):
    force = []
    for gap in gap_list:
        with open('gap.txt','w') as f:
            f.truncate()
            f.write('gap = '+ str(gap))
        print('APDL working... gap = '+ str(gap))
        os.system('apdl_3d.bat')
        with open('./output/force_y.txt','r') as f:
            force_y = float(f.readlines()[5][34:60])
            force.append(force_y)
        print('Done! Calculated Force = '+ str(force_y))
    return force
def analytic_sol(gap_list):
    force = []
    for gap in gap_list:
        C = eps_0*S/gap
        peak_coef = (2*pi*V_freq)*ind-1/(2*pi*V_freq)/C
        # print(V_freq*ind,1/V_freq/C) 
        q0 = V_amp/((2*pi*V_freq)*(res**2+(peak_coef)**2)**0.5)
        force.append( q0**2/(2*eps_0*S) )
    return force

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

def charge_ode(y, t):
     y1, y2 = y
     dydt = [y2, 1/ind*(V_amp*sin(2*pi*V_freq*t) - y1*init_gap/S/EPS_0 - res*y2)]
     return dydt

def calc_dynamic(init_gap,mass,time,dt):
    with open('params_pas_susp.txt','w') as f:
        f.truncate()
        f.write('mass = ' +    str(mass)+'\n')
        f.write('init_gap = '+ str(init_gap)  +'\n')
        f.write('time = ' +    str(time)+'\n')
        f.write('dt = ' +      str(dt)+'\n')
        f.write('NUM_CALCS = 1' +'\n')
        f.write("ANALYSIS_TYPE = 'dynamic_calc'")
    os.system('apdl_3d.bat')
    volt, time = read_from_file('./output/dyn_Volt_vs_Time.txt', True)
    return volt, time

def convergence(init_gap,mass,time,dt_list):
    v_conv = []
    for dt in dt_list:
        v_list, time_list = calc_dynamic(init_gap,mass,time,dt)
        v_conv.append(v_list[-1])
    return v_conv

def plot_convergence(dt_list, uy_conv):
    fig = plt.figure()
    fig.set_size_inches(13.5, 7)
    p1 = fig.add_subplot(111)
    p1.plot(dt_list,uy_conv,marker='o',linewidth=1, color='k',markersize=6, mew=6, linestyle='--')
    p1.set_xlabel("$Количество\ шагов\ интегрирования\ на\ период\ T=1/\omega$", fontsize=24)
    p1.set_ylabel("$Напряжение\ на\ конденасторе,\ В$", fontsize=24)
    p1.grid()
    fig.savefig('conv_plot',dpi=300)

# gap = 0.00048
# C = eps_0*S/gap
# phi0 = arctan(res/(V_freq*ind-1/V_freq/C))
# q0 = V_amp/((2*pi*V_freq)*(res**2+((2*pi*V_freq)*ind-1/(2*pi*V_freq)/C)**2)**0.5)

# q0_2 = V_amp*S*eps_0*4*pi/(  (ind*S*eps_0*4*pi*(2*pi*V_freq)**2-4*pi*gap)**2 + (res*S*eps_0*4*pi*(2*pi*V_freq))**2  )**0.5

# print('Capacitance:',C,'| R:',res,'| L:',ind,'| w:',V_freq)
# print('wL-1/wc',V_freq*ind -1/V_freq/C)
# print('phi0 =',phi0,'| sin(phi0)=',sin(phi0))
# print('q0 =',q0,'q0 =',q0_2)
# print('U_ic =',q0/C,'U_ic =',q0_2/C)
# print('xi =',1/2*res*(C/ind)**0.5)



# gap_list = list(arange(20e-6,500e-6,20e-6))
# # print(gap_list)


# force_an = analytic_sol(gap_list)
# # force_fe = get_force(gap_list)
# force_fe = [2.265521289e-06, 2.335369071e-06, 2.342284727e-06, 2.295555221e-06, 2.365718503e-06, 2.360149379e-06, 2.381396414e-06, 2.380296514e-06, 2.419122627e-06, 2.421736757e-06, 2.440887449e-06, 2.473360046e-06, 2.47786575e-06, 2.507216027e-06, 2.520965781e-06, 2.535945745e-06, 2.556265786e-06, 2.582918308e-06, 2.597908391e-06, 2.62077923e-06, 2.637391752e-06, 2.658381637e-06, 2.681000129e-06, 2.700353236e-06]
# print(force_fe)
# print(force_an)

# # compare = []
# # for i,j in zip(force_an,force_fe):
# #     k = i/j
# #     compare.append(k)
# # print(compare)


# fig = plt.figure()
# fig.set_size_inches(12, 12)
# p1 = fig.add_subplot(211)
# p1.plot(gap_list,force_fe,marker='o', label="",linewidth=2, color='r',markersize=6, mew=0)
# p2 = fig.add_subplot(212)
# p2.plot(gap_list,force_an,marker='o', label="",linewidth=1,linestyle='--', color='b',markersize=4, mew=0)
# filename = 'plot_force'
# fig.savefig(filename)

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)


# dt_list = [1/V_freq/(x*8) for x in range(1,12)]
# v_conv = convergence(init_gap,mass,1/V_freq*8,dt_list)
# v_conv = [-4.94092, -4.87135, -4.85313, -4.84711, -4.84467, -4.84355, -4.84302, -4.84277, -4.84266, -4.84262, -4.84263]
# dt_list = [1/V_freq/x for x in dt_list]
# print(v_conv)
# plot_convergence(dt_list, v_conv)

# for i in range(1,10):
#     print(1-v_conv[i]/v_conv[i-1])

# calc_dynamic(init_gap,mass,time,dt)

volt, time = read_from_file('./output/dyn_Volt_vs_Time.txt', True)
uy, time = read_from_file('./output/dyn_UY_vs_Time.txt', True)
force, time = read_from_file('./output/dyn_Force_vs_Time.txt', True)

# force = [x*253*2 for x in force]

C = EPS_0*S/init_gap
e_ode = odeint(charge_ode, [0.0, 0.0], time)[:, 0]
U_ode = [x/C for x in e_ode]
# U_ode = [x/C*0.863893771539 for x in e_ode]
F_ode = [-x**2/(2*EPS_0*S) for x in e_ode]

print(force[-1]/F_ode[-1])
print(volt[-1]/U_ode[-1]*100)

fig = plt.figure()
fig.set_size_inches(13.5, 21)
p1 = fig.add_subplot(311)
# p1.plot(time,uy,marker='o',linewidth=1, color='k',markersize=0, mew=0, linestyle='-')
# p1.set_xlabel("$Время,\ с$", fontsize=24)
# p1.set_ylabel("$Перемещение\ тела,\ м$", fontsize=24)
# p1.ticklabel_format(axis='both', style='sci', scilimits=(-2,2))
# p1.legend(loc='upper right',fontsize=22,numpoints=1)
# p1.grid()

p2 = fig.add_subplot(312)
p2.plot(time,volt,label='$FEM_{SOLID226}$',marker='o',linewidth=1, color='r',markersize=0, mew=0, linestyle='-')
p2.plot(time,U_ode,label='$ODE$',marker='o',linewidth=1, color='b',markersize=0, mew=0, linestyle='--')
p2.set_xlabel("$Время,\ с$", fontsize=24)
p2.set_ylabel("$Напряжение\ на\ электроде,\ В$", fontsize=24)
p2.ticklabel_format(axis='both', style='sci', scilimits=(-2,2))
p2.legend(loc='upper right',fontsize=22,numpoints=1)
p2.grid()

p3 = fig.add_subplot(313)
p3.plot(time,force,label='$FEM_{SOLID226}$',marker='o',linewidth=1, color='r',markersize=0, mew=0, linestyle='-')
p3.plot(time,F_ode,label='$ODE$',marker='o',linewidth=1, color='b',markersize=0, mew=0, linestyle='--')
p3.set_xlabel("$Время,\ с$", fontsize=24)
p3.set_ylabel("$Электрическая\ сила,\ Н$", fontsize=24)
p3.ticklabel_format(axis='both', style='sci', scilimits=(-2,2))
p3.legend(loc='upper right',fontsize=22,numpoints=1)
p3.grid()

fig.savefig('plot',dpi=300)    



# force_fe = [1.29338686, 0.323165135, 0.143615579, 0.08078252739, 0.05170096886, 0.0359037207, 0.02637846208, 0.02019616882, 0.01595758044]
