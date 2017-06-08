from numpy import pi, sin, sinh, cos, arange, subtract, around, exp, insert, arctan
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from config import *
import os
  
def read_from_file(fname, header = True, get_time_vector = False):
    '''
    Reads the results of capacitance calculation
    '''
    time_list, data_list = [], []
    with open(fname,'r') as f:
        if header == True:
            f = f.readlines()[14:]
        for line in f:
            time_list.append(float(line.split()[0]))
            data_list.append(float(line.split()[1]))
    if get_time_vector == True:
        return time_list, data_list
    else:
        return data_list

def ode_dynamic():

    time, fx_pos = read_from_file('./output/Fx_pos.txt')

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

def calc_cap_vs_gap(start_gap=1e-6, end_gap=20e-6, numcalc=20):
    '''
    Writes parameters of gap to file params.txt,
    starts ANSYS cap_calc.bat simulation and reads the results
    '''
    with open('cap_calc_params.txt','w') as f:
        f.truncate()
        f.write('d_min = '+     str(start_gap)+'\n')
        f.write('d_max = '+     str(end_gap)  +'\n')
        f.write('NUM_CALCS = '+ str(numcalc)  +'\n')
        f.write('d_inc = (d_max-d_min)/(NUM_CALCS-1)')

    os.system('capacitance_calc.bat')

    gap_list, cap_list = read_from_file('./output/CAPACITANCE_VS_GAP.txt', False, True)
    with open('cap_vs_gap_data.txt','w') as f:
        f.truncate()
        n = 0
        for cap, gap in zip(cap_list, gap_list):
            n += 1
            f.write('CAP'+str(n)+' = '+str(cap)+'\n')
            f.write('GAP'+str(n)+' = '+str(gap)+'\n')
    return cap_list, gap_list

def plot_cap():
    gap_list, cap_cae = read_from_file('./output/CAPACITANCE_VS_GAP.txt', False, True)
    cap_plate = plate_capacitance(gap_list)
    fig = plt.figure()
    fig.set_size_inches(10, 10)
    p1 = fig.add_subplot(111)
    p1.plot(gap_list,cap_cae,marker='o',label='$FEM$',linewidth=2, color='r',markersize=6, mew=0)
    p1.plot(gap_list,cap_plate,marker='o',label='$Plate$',linewidth=2, color='b',markersize=0, mew=0)
    p1.set_xlabel("$gap, m$", fontsize=15)
    p1.set_ylabel("$Capacinatce, F$", fontsize=15)
    p1.legend(loc='upper right',fontsize=15,numpoints=1)
    p1.grid()
    p1.set_title('$Capacitance$ $vs$ $gap$',fontsize=20)
    fig.savefig('capacitance_vs_gap',dpi=300)

def plate_capacitance(gap_list):
    cap_plate = []
    for gap in gap_list:
        cap_plate.append(eps0*S/gap)
    return cap_plate

def force_theory(ux_list, volt_list, volt_elec):
    force = []
    for x, v, v1 in zip(ux_list, volt_list, volt_elec):
        f = v**2*plate_capacitance([x+d])[0]/(x+d)/2 - v1**2*plate_capacitance([x-d])[0]/(x-d)/2
        force.append(f)
    return force


def main():
    '''Capacity versus gap calculation'''
    # calc_cap_vs_gap(10e-6,100e-6)
    # plot_cap()

    '''Launch ANSYS simulation'''
    start_simulation('suspension_3D.bat')

    '''Defining vectros to write data to'''
    u = [] 
    volt = []
    volt_in = []
    volt_rotor = []
    force = []
    sides = ["xp","xn", "yp","yn", "zp","zn"]

    '''Reading data from files'''
    time, bulk = read_from_file('./output/volt_rot.txt', True, True)
    u.append(read_from_file('./output/u_x.txt'))
    u.append(read_from_file('./output/u_y.txt'))
    u.append(read_from_file('./output/u_z.txt'))
    volt_in.append(read_from_file('./output/volt_xin.txt'))
    volt_in.append(read_from_file('./output/volt_yin.txt'))
    volt_in.append(read_from_file('./output/volt_zin.txt'))
    volt_rotor = read_from_file('./output/volt_rot.txt')

    for side in sides:
        volt.append(read_from_file('./output/volt_'+side+'.txt'))
        force.append(read_from_file('./output/force_'+side+'.txt'))

    volt_xsum, volt_ysum, volt_zsum = [], [], []
    force_xsum, force_ysum, force_zsum, = [], [], []
    for xp,xn,yp,yn,zp,zn in zip(force[0], force[1], force[2], force[3], force[4], force[5]):
        force_xsum.append(xp+xn)
        force_ysum.append(yp+yn)
        force_zsum.append(zp+zn)
    for vr,vx,vy,vz in zip(volt_rotor, volt[0], volt[2], volt[4]):
        volt_xsum.append(vr-vx)
        volt_ysum.append(vr-vy)
        volt_zsum.append(vr-vz)

    '''Plotting'''
    fig = plt.figure()
    fig.set_size_inches(12, 16)

    p1 = fig.add_subplot(411)
    p1.plot(time,u[0],label='$ux$',linewidth=1, color='r')
    p1.plot(time,u[1],label='$uy$',linewidth=1, color='g')
    p1.plot(time,u[2],label='$uz$',linewidth=1, color='b')
    p1.set_xlabel("$time, s$", fontsize=15)
    p1.set_ylabel("$u, m$", fontsize=15)
    p1.legend(loc='best',fontsize=15,numpoints=1)
    p1.grid()

    p2 = fig.add_subplot(412)
    p2.plot(time,volt_in[0],label='$volt-in-X$',linewidth=0.5, color='r')
    p2.plot(time,volt_in[1],label='$volt-in-Y$',linewidth=0.5, color='g')
    p2.plot(time,volt_in[2],label='$volt-in-Z$',linewidth=0.5, color='b')
    p2.plot(time,volt_rotor,label='$volt-rotor$', linewidth=1, color='y', linestyle='--')
    p2.plot(time,volt[0],label='$volt-electrode$',linewidth=1, color='c')
    p2.set_xlabel("$time, s$", fontsize=15)
    p2.set_ylabel("$volt, V$", fontsize=15)
    p2.legend(loc='upper right',fontsize=15,numpoints=1)
    p2.grid()

    p3 = fig.add_subplot(413)
    p3.plot(time,volt_zsum,label='$CapVolt-z$',linewidth=1, color='b')
    p3.plot(time,volt_xsum,label='$CapVolt-x$',linewidth=1, color='r')
    p3.plot(time,volt_ysum,label='$CapVolt-y$',linewidth=1, color='g')
    p3.set_xlabel("$time, s$", fontsize=15)
    p3.set_ylabel("$volt, V$", fontsize=15)
    p3.legend(loc='upper right',fontsize=15,numpoints=1)
    p3.grid()

    p4 = fig.add_subplot(414)
    p4.plot(time,force_zsum,label='$Zsum$',linewidth=2, color='b')
    p4.plot(time,force_xsum,label='$Xsum$',linewidth=2, color='r')
    p4.plot(time,force_ysum,label='$Ysum$',linewidth=2, color='g')
    p4.plot([time[0],time[-1]],[mass*g, mass*g],label='$mg$',linewidth=1, color='k')
    p4.set_xlabel("$time, s$", fontsize=15)
    p4.set_ylabel("$force, N$", fontsize=15)
    p4.legend(loc='upper right',fontsize=15,numpoints=1)
    p4.grid()

    fig.savefig('plot',dpi=300)

    # # fig2 = plt.figure()
    # # fig2.set_size_inches(12, 14)
    # # p1 = fig2.add_subplot(311)
    # # p1.plot(t_ux,ux,label='$ux$',linewidth=1, color='r')
    # # p1.set_xlabel("$time, s$", fontsize=15)
    # # p1.set_ylabel("$uy, m$", fontsize=15)
    # # p1.legend(loc='upper right',fontsize=15,numpoints=1)
    # # p1.grid()
    # # p2 = fig2.add_subplot(312)
    # # p2.plot(t_ux,uy,label='$uy$',linewidth=1, color='b')
    # # p2.set_xlabel("$time, s$", fontsize=15)
    # # p2.set_ylabel("$uy, m$", fontsize=15)
    # # p2.legend(loc='upper right',fontsize=15,numpoints=1)
    # # p2.grid()
    # # p3 = fig2.add_subplot(313)
    # # p3.plot(t_ux,uz,label='$uz$',linewidth=1, color='g')
    # # p3.set_xlabel("$time, s$", fontsize=15)
    # # p3.set_ylabel("$uz, m$", fontsize=15)
    # # p3.legend(loc='upper right',fontsize=15,numpoints=1)
    # # p3.grid()
    # # fig2.savefig('plot_u',dpi=300)

    # L = 0.08
    # R = 3636
    # w = 500e3
    # C = 5.5e-11
    # d0 = 40e-6

    # print(L*w, 1/C/w)
    # print(L*w - 1/C/w)
    # print(R/2/L)


if __name__ == "__main__":
    main()