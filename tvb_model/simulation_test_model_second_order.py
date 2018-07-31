import tvb.simulator.lab as lab
import numpy as np
import matplotlib.pylab as plt
import numpy.random as rgn
from Zerlaut_model_second_order import Zerlaut_model_second_order
# from Zerlaut import Zerlaut_second_order as Zerlaut_model_second_order
rgn.seed(10)
model = Zerlaut_model_second_order(variables_of_interest='E I C_ee C_ei C_ii'.split())
white_matter = lab.connectivity.Connectivity(load_default=True)
white_matter.speed = np.array([4.0])
white_matter_coupling = lab.coupling.Linear(a=0.0154)
heunint = lab.integrators.HeunDeterministic(dt=2**-6)
#Initialise some Monitors with period in physical time
mon_raw = lab.monitors.Raw(variables_of_interest='E I C_ee C_ei C_ii'.split())
mon_tavg = lab.monitors.TemporalAverage(variables_of_interest=[0,1,2,3,4],period=2**-2)

#Bundle them
what_to_watch = (mon_raw, mon_tavg)


#Initialise a Simulator -- Model, Connectivity, Integrator, and Monitors.
sim = lab.simulator.Simulator(model = model, connectivity = white_matter,
                          coupling = white_matter_coupling,
                          integrator = heunint, monitors = what_to_watch)

sim.configure()
# list = sim.run(simulation_length=2 ** 10)
# oscilator.P = np.ones((76))
# oscilator.P[0:5] = np.ones((5))*1000
# Perform the simulation
raw_data = []
raw_time = []
tavg_data = []
tavg_time = []

for raw, tavg in sim(simulation_length=2 ** 8):
    print("time",raw[0],"raw : ",np.mean(raw[1][0]),np.mean(raw[1][1]),np.mean(raw[1][2]),np.mean(raw[1][3]),np.mean(raw[1][4]))
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