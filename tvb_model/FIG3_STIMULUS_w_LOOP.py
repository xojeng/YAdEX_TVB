import tools_simulation
import parameter_76connectome_second_order_sleep
import numpy as np
import utils
import matplotlib.pyplot as plt

t=2000.0
parameter_76connectome_second_order_sleep.parameter_model['b_e'] = 100
parameter_76connectome_second_order_sleep.parameter_model['E_L_e'] = -55.0
parameter_76connectome_second_order_sleep.parameter_integrator['noise_parameter']['nsig']=[5e-9,5e-9,0.0,0.0,0.0,0.0,0.0,]
parameter_76connectome_second_order_sleep.parameter_coupling['parameter']['a']=0.02

Nregions = 76
parameter_stimulation = {}
parameter_stimulation['onset'] = 1500 # onset time in ms
parameter_stimulation['tau'] = 10 # stimulus duration? in ms
parameter_stimulation['T'] = 100000 # interstimulus interval in ms
parameter_stimulation['weights'] = np.zeros(Nregions)
parameter_stimulation['weights'][25] = 1e-1
#rPMCDL 25
#rM1 12
#lPMCDL 63

print(parameter_76connectome_second_order_sleep.parameter_model['b_e'], '\t', parameter_stimulation['weights'][25])
for seedy in range(1):#100):#100):

    name = 'adaptation'+str(parameter_76connectome_second_order_sleep.parameter_model['b_e'])+'_seed'+str(seedy)+'stim'+str()+'testplot'
    parameter_76connectome_second_order_sleep.parameter_simulation['path_result']= '/media/gpi/Elements/TVB_TA_neighbours_tstim1500/'+name+str(parameter_stimulation['weights'][25])


    #[0.7e-08, 0.7e-08, 0.0, 0.0, 0.0, 0.0, 0.0]

    simulator = tools_simulation.init(parameter_76connectome_second_order_sleep.parameter_simulation,
                                      parameter_76connectome_second_order_sleep.parameter_model,
                                      parameter_76connectome_second_order_sleep.parameter_connection_between_region,
                                      parameter_76connectome_second_order_sleep.parameter_coupling,
                                      parameter_76connectome_second_order_sleep.parameter_integrator,
                                      parameter_76connectome_second_order_sleep.parameter_monitor,
                                      parameter_stimulation = parameter_stimulation,
                                      my_seed = seedy)

    tools_simulation.run_simulation(simulator,
                                    t,
                                    parameter_76connectome_second_order_sleep.parameter_simulation,
                                    parameter_76connectome_second_order_sleep.parameter_monitor)

# tools_simulation.print_all(parameter_76connectome_second_order_sleep.parameter_simulation['path_result'],
#                            0.0,t,
#                            0,0)
#
# tools_simulation.print_all(parameter_76connectome_second_order_sleep.parameter_simulation['path_result'],
#                            0.0,t,
#                            0,5)
#
# tools_simulation.print_EI_one(parameter_76connectome_second_order_sleep.parameter_simulation['path_result'],
#                            0.0,t,
#                            0,25)

# # ##Plot into png for movies
# the_data = tools_simulation.get_result(parameter_76connectome_second_order_sleep.parameter_simulation['path_result'], 1000.0, t)
#
# for time_step in range(int(t*10)):
#     data  = the_data[0][1][time_step,0,utils.cortex.region_mapping]
#     utils.multiview(data, zmin=0.0, zmax = 0.02, shaded=False)
#     # plt.show()
#     plt.savefig(parameter_76connectome_second_order_sleep.parameter_simulation['path_result']+'/time_step_'+str(time_step)+'steps.png')
#     plt.close('all')


# plotting data average into single png
# the_data = tools_simulation.get_result(parameter_76connectome_second_order_sleep.parameter_simulation['path_result'], 1000.0, t)
# print(the_data[0][1].shape, t, utils.cortex.region_mapping)
# data = np.mean(the_data[0][1], axis = 0)[0,utils.cortex.region_mapping]
# print(data.shape)
# utils.multiview(data, zmin=0.0, zmax = 0.02, shaded=False)
#
# plt.show()
#
# plt.savefig(parameter_76connectome_second_order_sleep.parameter_simulation['path_result']+'/time_step_time_average_0to10Hzbar.png')
# # plt.close('all')