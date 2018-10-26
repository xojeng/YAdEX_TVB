import tvb.simulator.lab as lab
import numpy as np
import matplotlib.pylab as plt
import numpy.random as rgn
from Zerlaut import Zerlaut_adaptation_second_order as Zerlaut_model_second_order
import noise

rgn.seed(10)
# Initialisation of the model
model = Zerlaut_model_second_order(variables_of_interest='E I C_ee C_ei C_ii W_e W_i'.split())
model.state_variable_range['E'] = np.array([0.00,0.00])
model.state_variable_range['I'] = np.array([0.0,0.0])
model.state_variable_range['W_e'] = np.array([0.0,0.0])
model.state_variable_range['W_i'] = np.array([0.0,0.0])
model.E_L_e = -63.0
model.E_L_i = -65.0
# model.T = 10 # not need to change for the moment
model.external_input_ex_ex = 0.0
model.external_input_ex_in = 0.0
model.external_input_in_ex = 0.0
model.external_input_in_in = 0.0
# Configure the connectivity
white_matter = lab.connectivity.Connectivity(load_default=True)
white_matter.speed = np.array([1.0])
white_matter_coupling = lab.coupling.Difference(a=0.01)
# Deterministic integrator
##heunint = lab.integrators.HeunDeterministic(dt=5e-2)
# Stochastic integrator
hiss    = noise.Ornstein_Ulhenbeck_process(T=10.0,nsig=np.array([0.08e-4,0.5e-4,0.0,0.0,0.0,0.0,0.0]), ntau=1.0)
heunint = lab.integrators.HeunStochastic(dt=5e-2, noise=hiss)
#Initialise some Monitors with period in physical time ms and interes variable
mon_raw = lab.monitors.Raw(variables_of_interest='E I C_ee C_ei C_ii W_e W_i'.split())    # value of the model
mon_tavg = lab.monitors.TemporalAverage(variables_of_interest=[0,1,2,3,4,5,6],period=1.0) # average in time of value
# mon_bold = lab.monitors.Bold(variables_of_interest=[0],period=2000)         #in ms => simulated bold signal

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

for raw, tavg in sim(simulation_length=20000.0):
    if not raw is None:
        print("time",raw[0],"max/mean : ",np.max(raw[1][0]),np.mean(raw[1][0]),np.max(raw[1][1]),np.mean(raw[1][1]),np.max(raw[1][5]),np.mean(raw[1][5]))
        raw_time.append(raw[0])
        raw_data.append(raw[1])

    if not tavg is None:
        tavg_time.append(tavg[0])
        tavg_data.append(tavg[1])

#Make the lists numpy.arrays for easier use.
raw_data = np.array(raw_data)
np.save('./RAW.npy',raw_data)
tavg_data = np.array(tavg_data)
np.save('./TAVG.npy',tavg_data)

# #Plot raw time series
# plt.figure(1)
# plt.plot(raw_time, raw_data[:, 0, :, 0])
# plt.title("Raw -- State variable E")
# plt.show()
# plt.close('all')
# plt.figure(2)
# plt.plot(raw_time, raw_data[:, 1, :, 0])
# plt.title("Raw -- State variable I")
# plt.show()
# plt.close('all')
# plt.figure(3)
# plt.plot(raw_time, raw_data[:, 2, :, 0])
# plt.title("Raw -- State variable C_EE")
# plt.show()
# plt.close('all')
# plt.figure(4)
# plt.plot(raw_time, raw_data[:, 3, :, 0])
# plt.title("Raw -- State variable C_EI")
# plt.show()
# plt.close('all')
# plt.figure(5)
# plt.plot(raw_time, raw_data[:, 4, :, 0])
# plt.title("Raw -- State variable C_II")
# plt.show()
# plt.close('all')
# plt.figure(6)
# plt.plot(raw_time, raw_data[:, 5, :, 0])
# plt.title("Raw -- State variable W_e")
# plt.show()
# plt.close('all')
# plt.figure(7)
# plt.plot(raw_time, raw_data[:, 6, :, 0])
# plt.title("Raw -- State variable W_i")
# plt.show()
# plt.close('all')

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
# plt.show()
