from numpy import pi, sin, sinh, cos, arange, subtract, around, exp, insert,arctan
import matplotlib.pyplot as plt
from config import *
import os

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

gap = 0.00048
C = eps_0*S/gap
phi0 = arctan(res/(V_freq*ind-1/V_freq/C))
q0 = V_amp/((2*pi*V_freq)*(res**2+((2*pi*V_freq)*ind-1/(2*pi*V_freq)/C)**2)**0.5)

q0_2 = V_amp*S*eps_0*4*pi/(  (ind*S*eps_0*4*pi*(2*pi*V_freq)**2-4*pi*gap)**2 + (res*S*eps_0*4*pi*(2*pi*V_freq))**2  )**0.5

print('Capacitance:',C,'| R:',res,'| L:',ind,'| w:',V_freq)
print('wL-1/wc',V_freq*ind -1/V_freq/C)
print('phi0 =',phi0,'| sin(phi0)=',sin(phi0))
print('q0 =',q0,'q0 =',q0_2)
print('U_ic =',q0/C,'U_ic =',q0_2/C)
print('xi =',1/2*res*(C/ind)**0.5)



gap_list = list(arange(20e-6,500e-6,20e-6))
print(gap_list)


force_an = analytic_sol(gap_list)
# force_fe = get_force(gap_list)
force_fe = [2.265521289e-06, 2.335369071e-06, 2.342284727e-06, 2.295555221e-06, 2.365718503e-06, 2.360149379e-06, 2.381396414e-06, 2.380296514e-06, 2.419122627e-06, 2.421736757e-06, 2.440887449e-06, 2.473360046e-06, 2.47786575e-06, 2.507216027e-06, 2.520965781e-06, 2.535945745e-06, 2.556265786e-06, 2.582918308e-06, 2.597908391e-06, 2.62077923e-06, 2.637391752e-06, 2.658381637e-06, 2.681000129e-06, 2.700353236e-06]
print(force_fe)
print(force_an)

# # compare = []
# # for i,j in zip(force_an,force_fe):
# #     k = i/j
# #     compare.append(k)
# # print(compare)

fig = plt.figure()
fig.set_size_inches(12, 12)
p1 = fig.add_subplot(211)
p1.plot(gap_list,force_fe,marker='o', label="",linewidth=2, color='r',markersize=6, mew=0)
p2 = fig.add_subplot(212)
p2.plot(gap_list,force_an,marker='o', label="",linewidth=1,linestyle='--', color='b',markersize=4, mew=0)
filename = 'plot_force'
fig.savefig(filename)





# force_fe = [1.29338686, 0.323165135, 0.143615579, 0.08078252739, 0.05170096886, 0.0359037207, 0.02637846208, 0.02019616882, 0.01595758044]
