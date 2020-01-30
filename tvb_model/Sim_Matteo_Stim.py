import tools_simulation
import parameter_Matteo
import numpy as np
import utils
import matplotlib.pyplot as plt

t=14000.0

parameter_Matteo.parameter_simulation['path_result']= '/mnt/usb-WD_easystore_25FC_575856314136373548355350-0:0-part1/4TB_copy/Athens_Matteo_TVB_params/MTOO_b0_try/'
    #'/media/gpi/4TB/Athens_Matteo_TVB_params/MTOO_b100_try'

Nregions = 84
parameter_stimulation = {}
parameter_stimulation['onset'] = 9000 # onset time in ms
parameter_stimulation['tau'] = 10 # stimulus duration in ms
parameter_stimulation['T'] = 100000 # interstimulus interval in ms
parameter_stimulation['weights'] = np.zeros(Nregions)
parameter_stimulation['weights'][59] = 1e-3

simulator = tools_simulation.init(parameter_Matteo.parameter_simulation,
                                 parameter_Matteo.parameter_model,
                                 parameter_Matteo.parameter_connection_between_region,
                                 parameter_Matteo.parameter_coupling,
                                 parameter_Matteo.parameter_integrator,
                                 parameter_Matteo.parameter_monitor,
                                 parameter_stimulation = parameter_stimulation)

# tools_simulation.run_simulation(simulator,
#                                t,
#                                parameter_Matteo.parameter_simulation,
#                                parameter_Matteo.parameter_monitor)

tools_simulation.print_all(parameter_Matteo.parameter_simulation['path_result'],
                           4000.0,t,
                           0,simulator.connectivity)

# tools_simulation.print_all(parameter_Matteo.parameter_simulation['path_result'],
#                            1.0,t,
#                            0,simulator.connectivity)

tools_simulation.print_EI_one(parameter_Matteo.parameter_simulation['path_result'],
                             4000.0,t,
                             0,59)

# the_data = tools_simulation.get_result(parameter_Matteo.parameter_simulation['path_result'], 0.0, t)
#
# for time_step in range(int(t*10)):
#     data  = the_data[0][1][time_step,0,utils.cortex.region_mapping]
#     utils.multiview(data, zmin=0.0, zmax = 0.02, shaded=False)
#     # plt.show()
#     plt.savefig(parameter_Matteo.parameter_simulation['path_result']+'/time_step_'+str(time_step)+'.png')
#     plt.close('all')

