import logging
logging.getLogger('numba').setLevel('ERROR')
logging.getLogger('matplotlib').setLevel('ERROR')
logging.getLogger('tvb').setLevel('ERROR')
import tools_simulation
import  parameter_New_1 as parameter
import numpy as np
import utils
import matplotlib.pyplot as plt

t=200.0

parameter.parameter_simulation['path_result']='./new_parameter_test_stim'

# tools_simulation.print_bistability(parameter.parameter_model)

simulator = tools_simulation.init(parameter.parameter_simulation,
                                 parameter.parameter_model,
                                 parameter.parameter_connection_between_region,
                                 parameter.parameter_coupling,
                                 parameter.parameter_integrator,
                                 parameter.parameter_monitor,
                                  parameter_stimulation={'onset':100.0,"tau":10.0,"T":100.0,"weights":np.ones(84)})
tools_simulation.run_simulation(simulator,
                               t,
                               parameter.parameter_simulation,
                               parameter.parameter_monitor)

tools_simulation.print_all(parameter.parameter_simulation['path_result'],
                          0.0,t,
                          0,simulator.connectivity)
#
# tools_simulation.print_all(parameter_PosterWeek_desperation_disconnected.parameter_simulation['path_result'],
#                           0.0,t,
#                           0,5)

tools_simulation.print_EI_one(parameter.parameter_simulation['path_result'],
                             0.0,t,
                            0,0)
tools_simulation.print_EI_one(parameter.parameter_simulation['path_result'],
                             0.0,t,
                            0,21)
# the_data = tools_simulation.get_result(parameter.parameter_simulation['path_result'], 0.0, t)
# utils.multiview(the_data[0][1][0,0,utils.cortex.region_mapping], zmin=0.0, zmax = 0.02, shaded=False)

# for time_step in range(int(t*10)):
#     data  = the_data[0][1][time_step,0,utils.cortex.region_mapping]
#     utils.multiview(data, zmin=0.0, zmax = 0.02, shaded=False)
#     # plt.show()
#     plt.savefig(parameter_PosterWeek_desperation_disconnected.parameter_simulation['path_result']+'/time_step_'+str(time_step)+'.png')
#     plt.close('all')

