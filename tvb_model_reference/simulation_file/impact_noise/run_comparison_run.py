from tvb_model_reference.simulation_file.impact_noise.parameter_TVB import Parameter
import tvb_model_reference.src.tools_simulation as tools
import matplotlib.pyplot as plt
import copy

parameters = Parameter()
# tools.print_bistability(parameters.parameter_model,show=False)
run_sim = 10000.0
simulator = tools.init(parameters.parameter_simulation,
                              parameters.parameter_model,
                              parameters.parameter_connection_between_region,
                              parameters.parameter_coupling,
                              parameters.parameter_integrator,
                              parameters.parameter_monitor)
tools.run_simulation(simulator,
                            run_sim,
                            parameters.parameter_simulation,
                            parameters.parameter_monitor)
result = tools.get_result(parameters.parameter_simulation['path_result'],0.0,run_sim)
ax1 = plt.subplot(221)
ax1.plot(result[0][0]*1e-3,result[0][1][:,0,0]*1e3)
ax1.plot(result[0][0]*1e-3,result[0][1][:,1,0]*1e3)
ax1.set_title('firing rate')
ax1.set_xlabel('time in s')
ax1.set_ylabel('Mean firing rate in Hz')
ax2 = plt.subplot(222)
ax2.plot(result[0][0]*1e-3,result[0][1][:,5,0])
ax2.set_title('adaptation rate')
ax2.set_xlabel('time in s')
ax2.set_ylabel('nA')
ax3 = plt.subplot(223)
ax3.plot(result[0][0]*1e-3,(result[0][1][:,7,0]*parameters.parameter_model['weight_noise']+parameters.parameter_model['external_input_ex_ex'])*1e3)
ax3.set_title('noise')
ax3.set_xlabel('time in s')
ax3.set_ylabel('Hz ')
plt.show()
