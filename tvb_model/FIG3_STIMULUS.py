import tools_simulation
import parameter_76connectome_second_order_sleep
import numpy as np
import utils
import matplotlib.pyplot as plt

#change the save path before running!
parameter_76connectome_second_order_sleep.parameter_simulation['path_result']= '/home/gpi/Desktop/TVB_Linux_1.5.4/Lionel_updates_June/YAdEX_TVB-master/tvb_model/Fig2_AI_60s/'

t=60000.0
#change b_e to transition between states. b_e=0 for AI-like state and b_e=100 for SWS-like state
parameter_76connectome_second_order_sleep.parameter_model['b_e'] = 100.0
parameter_76connectome_second_order_sleep.parameter_model['E_L_e'] = -55.0
parameter_76connectome_second_order_sleep.parameter_integrator['noise_parameter']['nsig']=[5e-9,5e-9,0.0,0.0,0.0,0.0,0.0,]
parameter_76connectome_second_order_sleep.parameter_coupling['parameter']['a']=0.02
#
Nregions = 76
parameter_stimulation = {}
parameter_stimulation['onset'] = 6000 # onset time in ms
parameter_stimulation['tau'] = 500 # stimulus duration? in ms
parameter_stimulation['T'] = 50000 # interstimulus interval in ms
parameter_stimulation['weights'] = np.zeros(Nregions)

#regions of stimulus:
#rPMCDL 25
#rM1 12
#lPMCDL 63
parameter_stimulation['weights'][25] = 1e-3

simulator = tools_simulation.init(parameter_76connectome_second_order_sleep.parameter_simulation,
                                  parameter_76connectome_second_order_sleep.parameter_model,
                                  parameter_76connectome_second_order_sleep.parameter_connection_between_region,
                                  parameter_76connectome_second_order_sleep.parameter_coupling,
                                  parameter_76connectome_second_order_sleep.parameter_integrator,
                                  parameter_76connectome_second_order_sleep.parameter_monitor,
                                  parameter_stimulation = parameter_stimulation,
                                  my_seed = 1)

tools_simulation.run_simulation(simulator,
                                t,
                                parameter_76connectome_second_order_sleep.parameter_simulation,
                                parameter_76connectome_second_order_sleep.parameter_monitor)

tools_simulation.print_all(parameter_76connectome_second_order_sleep.parameter_simulation['path_result'],
                           0.1000,t,
                           0,0)

tools_simulation.print_all(parameter_76connectome_second_order_sleep.parameter_simulation['path_result'],
                           0.1000,t,
                           0,0)

tools_simulation.print_EI_one(parameter_76connectome_second_order_sleep.parameter_simulation['path_result'],
                           0.1000,t,
                           0,25)

# # ##Plot into pngs; one png for each time step for movies
# the_data = tools_simulation.get_result(parameter_76connectome_second_order_sleep.parameter_simulation['path_result'], 1000.0, t)
#
# for time_step in range(int(t*10)):
#     data  = the_data[0][1][time_step,0,utils.cortex.region_mapping]
#     utils.multiview(data, zmin=0.0, zmax = 0.02, shaded=False)
#     # plt.show()
#     plt.savefig(parameter_76connectome_second_order_sleep.parameter_simulation['path_result']+'/time_step_'+str(time_step)+'steps.png')
#     plt.close('all')

the_data = tools_simulation.get_result(parameter_76connectome_second_order_sleep.parameter_simulation['path_result'], 1000.0, t)
data  = the_data[0][1][0,0,utils.cortex.region_mapping]
print(data.shape)
utils.multiview(data, zmin=0.0, zmax = 0.02, shaded=False)
plt.show()
#     plt.savefig(parameter_76connectome_second_order_sleep.parameter_simulation['path_result']+'/time_step_'+str(time_step)+'steps.png')
#     plt.close('all')

# # plotting data average into single png
the_data = tools_simulation.get_result(parameter_76connectome_second_order_sleep.parameter_simulation['path_result'], 1000.0, t)
data = np.mean(the_data[0][1], axis = 0)[0,utils.cortex.region_mapping]
print(data)
print(data.shape)
utils.multiview(data, zmin=0.0, zmax = 0.02, shaded=False)
plt.savefig(parameter_76connectome_second_order_sleep.parameter_simulation['path_result']+'/time_step_time_average_0to10Hzbar.png')
plt.close('all')