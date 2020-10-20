import tvb.simulator.lab as lab
import numpy as np
import matplotlib.pylab as plt
import numpy.random as rgn
from Zerlaut import Zerlaut_adaptation_second_order as Zerlaut_model_second_order
import noise



rgn.seed(10)
# Initialisation of the model
model = Zerlaut_model_second_order(variables_of_interest='E I C_ee C_ei C_ii W_e W_i'.split())
# init default 0
# model.state_variable_range['E'] = np.array([1.0e-1,0.5e-3])
# model.state_variable_range['I'] = np.array([1.0e-1,0.5e-3])
# model.state_variable_range['W_e'] = np.array([100.0,10.0])
# model.state_variable_range['W_i'] = np.array([0.0,0.0])

#init 2
# model.state_variable_range['E'] = np.array([1.0e-2,0.5e-3])
# model.state_variable_range['I'] = np.array([1.0e-2,0.5e-3])
# model.state_variable_range['W_e'] = np.array([100.0,10.0])
# model.state_variable_range['W_i'] = np.array([0.0,0.0])

#init 3
# model.state_variable_range['E'] = np.array([1.0e-2,0.5e-3])
# model.state_variable_range['I'] = np.array([1.0e-2,0.5e-3])
# model.state_variable_range['W_e'] = np.array([1000.0,500.0])
# model.state_variable_range['W_i'] = np.array([0.0,0.0])

#init 4
model.state_variable_range['E'] = np.array([0.0e-2,0e-3])
model.state_variable_range['I'] = np.array([0.0e-2,0e-3])
model.state_variable_range['W_e'] = np.array([1000.0,500.0])
model.state_variable_range['W_i'] = np.array([0.0,0.0])

#init 5
# model.state_variable_range['E'] = np.array([0.0e-2,0e-3])
# model.state_variable_range['I'] = np.array([0.0e-2,0e-3])
# model.state_variable_range['W_e'] = np.array([100.0,0.0])
# model.state_variable_range['W_i'] = np.array([0.0,0.0])

# Configure the connectivity
white_matter = lab.connectivity.Connectivity(load_default=False)
white_matter.speed = np.array([1.0])
white_matter_coupling = lab.coupling.Linear(a=0.05)
# Deterministic integrator
##heunint = lab.integrators.HeunDeterministic(dt=5e-2)
# Stochastic integrator
hiss    = noise.Ornstein_Ulhenbeck_process(nsig=np.array([5.0e-6,5.0e-6,0.0,0.0,0.0,0.0,0.0]))
heunint = lab.integrators.HeunStochastic(dt=5e-2, noise=hiss)
# heunint = lab.integrators.HeunDeterministic(dt=5e-2)
#Initialise some Monitors with period in physical time ms and interes variable
mon_raw = lab.monitors.Raw(variables_of_interest='E I C_ee C_ei C_ii W_e W_i'.split())    # value of the model
mon_tavg = lab.monitors.TemporalAverage(variables_of_interest=[0,1,2,3,4,5,6],period=5.0) # average in time of value
# mon_bold = lab.monitors.Bold(variables_of_interest=[0],period=2000)         #in ms => simulated bold signal

#Bundle them
what_to_watch = (mon_raw, mon_tavg)
# what_to_watch = (mon_bold)

#Initialise a Simulator -- Model, Connectivity, Integrator, and Monitors.
sim = lab.simulator.Simulator(model = model, connectivity = white_matter,
                          coupling = white_matter_coupling,
                          integrator = heunint, monitors = what_to_watch)

# sim.configure()
# bold_time =[]
# bold_data =[]
# for bold in sim(simulation_length=1200000.0):
#     if not bold is None:
#         print("time",bold[0][0],"max/mean : ",np.max(bold[0][1][0]),np.mean(bold[0][1][0]),np.min(bold[0][1][0]))
#         bold_time.append(bold[0][0])
#         bold_data.append(bold[0][1])
# Make the lists numpy.arrays for easier use.
# bold_data = np.array(bold_data)
# np.save('./test_20/BOLD_end.npy',bold_data)
# bold_time = np.array(bold_time)
# np.save('./test_20/BOLDT_end.npy',bold_time)



sim.configure()
raw_data = []
raw_time = []
tavg_data = []
tavg_time = []

for raw, tavg in sim(simulation_length=80000.0):
    if not raw is None:
        raw_time.append(raw[0])
        raw_data.append(raw[1])

    if not tavg is None:
        print("time",tavg[0],"max/mean : ",np.max(tavg[1][0]),np.mean(tavg[1][0]),np.max(tavg[1][1]),np.mean(tavg[1][1]),np.max(tavg[1][5]),np.mean(tavg[1][5]))
        tavg_time.append(tavg[0])
        tavg_data.append(tavg[1])

#Make the lists numpy.arrays for easier use.
raw_data = np.array(raw_data)
np.save('./test_good_noise/noise_0_05/RAW_5e_6_init_4.npy',raw_data)
tavg_data = np.array(tavg_data)
np.save('./test_good_noise/noise_0_05/TAVG_5e_6_init_4.npy',tavg_data)

#Plot raw time series

# plt.figure(1)
# plt.plot(raw_time, raw_data[:, 0, :, 0])
# plt.title("Raw -- State variable E")
# plt.show()
# plt.close('all')
# plt.figure(2)
# plt.plot(raw_time, raw_data[:, 1, :, 0])
# plt.title("Raw -- State variable I")
# # plt.show()
# # plt.close('all')
# plt.figure(3)
# plt.plot(raw_time, raw_data[:, 2, :, 0])
# plt.title("Raw -- State variable C_EE")
# # plt.show()
# # plt.close('all')
# plt.figure(4)
# plt.plot(raw_time, raw_data[:, 3, :, 0])
# plt.title("Raw -- State variable C_EI")
# # plt.show()
# # plt.close('all')
# plt.figure(5)
# plt.plot(raw_time, raw_data[:, 4, :, 0])
# plt.title("Raw -- State variable C_II")
# # plt.show()
# # plt.close('all')
# plt.figure(6)
# plt.plot(raw_time, raw_data[:, 5, :, 0])
# plt.title("Raw -- State variable W_e")
# # plt.show()
# # plt.close('all')
# plt.figure(7)
# plt.plot(raw_time, raw_data[:, 6, :, 0])
# plt.title("Raw -- State variable W_i")
# plt.show()
# # plt.close('all')

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
# # Show them
# plt.show()
