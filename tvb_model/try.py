import parameter_oneRegion_second_order_Wake_mouse as parameter
import tools_simulation
import numpy as np
import tools_plotting

rangeG = [0.05, 0.1]
path = parameter.parameter_simulation['path_result']+'high_G_0.1/'
parameter.parameter_simulation['path_result'] +="test/"

# tools_simulation.gen_bif_data(rangeG, parameter.parameter_simulation,
#                                   parameter.parameter_model,
#                                   parameter.parameter_connection_between_region,
#                                   parameter.parameter_coupling,
#                                   parameter.parameter_integrator,
#                                   parameter.parameter_monitor)

# result = tools_simulation.get_result(path, 4999., 5000.)[0][1][-1]
# print(result)

tools_plotting.plot_bif_diagram(rangeG, parameter.parameter_simulation, save=True, plot_name="test")