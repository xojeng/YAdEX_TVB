import tools_simulation
import parameter_PosterWeek_desperation_disconnected
import numpy as np
import utils
import matplotlib.pyplot as plt

t=10000.0

parameter_PosterWeek_desperation_disconnected.parameter_simulation['path_result']='./posterweek_apoint2_b0_disconnected'

simulator = tools_simulation.init(parameter_PosterWeek_desperation_disconnected.parameter_simulation,
                                 parameter_PosterWeek_desperation_disconnected.parameter_model,
                                 parameter_PosterWeek_desperation_disconnected.parameter_connection_between_region,
                                 parameter_PosterWeek_desperation_disconnected.parameter_coupling,
                                 parameter_PosterWeek_desperation_disconnected.parameter_integrator,
                                 parameter_PosterWeek_desperation_disconnected.parameter_monitor)
tools_simulation.run_simulation(simulator,
                               t,
                               parameter_PosterWeek_desperation_disconnected.parameter_simulation,
                               parameter_PosterWeek_desperation_disconnected.parameter_monitor)

# tools_simulation.print_all(parameter_PosterWeek_desperation_disconnected.parameter_simulation['path_result'],
#                           0.0,t,
#                           0,0)
#
# tools_simulation.print_all(parameter_PosterWeek_desperation_disconnected.parameter_simulation['path_result'],
#                           0.0,t,
#                           0,5)

tools_simulation.print_EI_one(parameter_PosterWeek_desperation_disconnected.parameter_simulation['path_result'],
                             0.0,t,
                            0,0)

# the_data = tools_simulation.get_result(parameter_PosterWeek_desperation_disconnected.parameter_simulation['path_result'], 0.0, t)
#
# for time_step in range(int(t*10)):
#     data  = the_data[0][1][time_step,0,utils.cortex.region_mapping]
#     utils.multiview(data, zmin=0.0, zmax = 0.02, shaded=False)
#     # plt.show()
#     plt.savefig(parameter_PosterWeek_desperation_disconnected.parameter_simulation['path_result']+'/time_step_'+str(time_step)+'.png')
#     plt.close('all')

