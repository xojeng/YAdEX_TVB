# import os; os.chdir('C:\Users\internet\Documents\doctora\\travail\YAdEX_TVB\\tvb_model')
import tvb.simulator.lab as lab
import numpy as np
import matplotlib.pylab as plt
import numpy.random as rgn
# from Zerlaut_model_second_order import Zerlaut_model_second_order
from Zerlaut import Zerlaut_adaptation_second_order as Zerlaut_model_second_order

rgn.seed(10)
# Initialisation of the model
model = Zerlaut_model_second_order(variables_of_interest='E I C_ee C_ei C_ii W_e W_i'.split())
##model.state_variable_range['E'] = np.array([1.4907357840560944e-3,1.4907357840560944e-3])
##model.state_variable_range['I'] = np.array([9.198083797700072e-3,9.198083797700072e-3])
##model.state_variable_range['C_EE'] = np.array([0.15129606503776358e-3,0.15129606503776358e-3])
##model.state_variable_range['C_EI'] = np.array([0.3025868361216362e-3,0.3025868361216362e-3])
##model.state_variable_range['C_II'] = np.array([0.7632335165277322e-3,0.7632335165277322e-3])
##model.state_variable_range['W_e'] = np.array([7.454102191637575e1,7.454102191637575e1])
##model.state_variable_range['W_i'] = np.array([0.0,0.0])
# model.state_variable_range['E'] = np.array([2.8e-3,2.8e-3])
# model.state_variable_range['I'] = np.array([50.0e-3,50.0e-3])
# model.state_variable_range['C_EE'] = np.array([0.5e-3,0.5e-3])
# model.state_variable_range['C_EI'] = np.array([0.5e-3,0.5e-3])
# model.state_variable_range['C_II'] = np.array([0.5e-3,0.5e-3])
# print(2.8e-3*model.b_e*model.tau_w_e)
# model.state_variable_range['W_e'] = np.array([2.8e-3*model.b_e*model.tau_w_e,2.8e-3*model.b_e*model.tau_w_e])
# model.state_variable_range['W_i'] = np.array([0.0,0.0])
# model.external_input_ex_ex = 2.5e-3
# model.external_input_ex_in = 0.0
# model.external_input_in_ex = 2.5e-3
# model.external_input_in_in = 0.0

# model.state_variable_range['E'] = np.array([0.0e-3,0.0e-3])
# model.state_variable_range['I'] = np.array([0.0e-3,0.0e-3])
# model.state_variable_range['C_EE'] = np.array([0.0e-3,0.0e-3])
# model.state_variable_range['C_EI'] = np.array([0.0e-3,0.0e-3])
# model.state_variable_range['C_II'] = np.array([0.0e-3,0.0e-3])
# model.state_variable_range['E'] = np.array([1.02128531393361e-003,1.02128531393361e-003])
# model.state_variable_range['I'] = np.array([1.84631123578572e-003,1.84631123578572e-003])
# model.state_variable_range['C_EE'] = np.array([3.12954495998902e-009,3.12954495998902e-009])
# model.state_variable_range['C_EI'] = np.array([-27.8312775414693e-012,-27.8312775414693e-012])
# model.state_variable_range['C_II'] = np.array([22.1680017413319e-009,22.1680017413319e-009])
model.a_e =0.0
model.b_e = 0.0
model.a_i =0.0
model.b_i =0.0
model.tau_w_e = 1.0
model.tau_w_i =1.0
model.E_L_e = -65.2
model.E_L_i = -65.2
model.state_variable_range['W_e'] = np.array([0.0e-3,0.0e-3])
model.state_variable_range['W_i'] = np.array([0.0,0.0])
model.P_e = [-49.8e-3,5.06e-3,-25.0e-3,1.4e-3,-0.41e-3,10.5e-3,-36.0e-3,7.4e-3,1.2e-3,-40.7e-3]
model.P_i = [-51.4e-3,4.0e-3,-8.3e-3,0.2e-3,-0.5e-3,1.4e-3,-14.6e-3,4.5e-3,2.8e-3,-15.3e-3]
# Configure the connectivity
white_matter = lab.connectivity.Connectivity(number_of_regions=2,tract_lengths=np.array([[1.0,1.0],[1.0,1.0]]),weights=np.array([[0.0,0.0],[0.0,0.0]]))
white_matter.speed = np.array([4.0])
white_matter_coupling = lab.coupling.Linear(a=0.0)
# Deterministic integrator
##heunint = lab.integrators.HeunDeterministic(dt=2**-6)
heunint = lab.integrators.HeunDeterministic(dt=5e-3)
# Stochastic integrator
# hiss    = lab.noise.Multiplicative(nsig = np.array([0.0015]))
# heunint = lab.integrators.HeunStochastic(dt=2**-6, noise=hiss)
#Initialise some Monitors with period in physical time
mon_raw = lab.monitors.Raw(variables_of_interest='E I C_ee C_ei C_ii W_e W_i'.split())
mon_tavg = lab.monitors.TemporalAverage(variables_of_interest=[0,1,2,3,4,5,6],period=2**-2)

#Bundle them
what_to_watch = (mon_raw, mon_tavg)

#Initialise a Simulator -- Model, Connectivity, Integrator, and Monitors.
sim = lab.simulator.Simulator(model = model, connectivity = white_matter,
                          coupling = white_matter_coupling,
                          integrator = heunint, monitors = what_to_watch)

sim.configure()
raw_data = []
raw_time = []
tavg_data = []
tavg_time = []

model.external_input_ex_ex = 0.0e-3
model.external_input_in_ex = 0.0e-3
for raw, tavg in sim(simulation_length=200.0):
    print("time",raw[0],"raw : ",np.mean(raw[1][0]),np.mean(raw[1][1]),np.mean(raw[1][2]),np.mean(raw[1][3]),np.mean(raw[1][4]),np.mean(raw[1][5]),np.mean(raw[1][6]))
    if not raw is None:
        raw_time.append(raw[0])
        raw_data.append(raw[1])

    if not tavg is None:
        tavg_time.append(tavg[0])
        tavg_data.append(tavg[1])

model.external_input_ex_ex = 0.0e-3
model.external_input_in_ex = 0.0e-3
for raw, tavg in sim(simulation_length=500.0):
    print("time",raw[0],"raw : ",np.mean(raw[1][0]),np.mean(raw[1][1]),np.mean(raw[1][2]),np.mean(raw[1][3]),np.mean(raw[1][4]),np.mean(raw[1][5]),np.mean(raw[1][6]))
    if not raw is None:
        raw_time.append(raw[0])
        raw_data.append(raw[1])

    if not tavg is None:
        tavg_time.append(tavg[0])
        tavg_data.append(tavg[1])

#Make the lists numpy.arrays for easier use.
RAW = np.array(raw_data)
TAVG = np.array(tavg_data)

#Plot raw time series
plt.figure(1)
plt.plot(raw_time, RAW[:, 0, :, 0])
plt.title("Raw -- State variable E")
plt.figure(2)
plt.plot(raw_time, RAW[:, 1, :, 0])
plt.title("Raw -- State variable I")
plt.figure(3)
plt.plot(raw_time, RAW[:, 2, :, 0])
plt.title("Raw -- State variable C_EE")
plt.figure(4)
plt.plot(raw_time, RAW[:, 3, :, 0])
plt.title("Raw -- State variable C_EI")
plt.figure(5)
plt.plot(raw_time, RAW[:, 4, :, 0])
plt.title("Raw -- State variable C_II")
plt.figure(6)
plt.plot(raw_time, RAW[:, 5, :, 0])
plt.title("Raw -- State variable W_e")
plt.figure(7)
plt.plot(raw_time, RAW[:, 6, :, 0])
plt.title("Raw -- State variable W_i")

# #Plot temporally averaged time series
# plt.figure(6)
# plt.plot(tavg_time, TAVG[:, 0, :, 0])
# plt.title("Temporal average-- State variable E")
# plt.figure(7)
# plt.plot(tavg_time, TAVG[:, 1, :, 0])
# plt.title("Temporal average-- State variable I")
# plt.figure(8)
# plt.plot(tavg_time, TAVG[:, 2, :, 0])
# plt.title("Temporal average-- State variable C_EE")
# plt.figure(9)
# plt.plot(tavg_time, TAVG[:, 3, :, 0])
# plt.title("Temporal average-- State variable C_EI")
# plt.figure(10)
# plt.plot(tavg_time, TAVG[:, 4, :, 0])
# plt.title("Temporal average-- State variable C_II")
#Show them
plt.show()
