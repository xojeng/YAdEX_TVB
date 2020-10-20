import tools_simulation
import parameter_76connectome_second_order_sleep
import numpy as np

parameter_76connectome_second_order_sleep.parameter_simulation['path_result']='./76connectome_secondorder_sleep3/'

simulator = tools_simulation.init(parameter_76connectome_second_order_sleep.parameter_simulation,
                                  parameter_76connectome_second_order_sleep.parameter_model,
                                  parameter_76connectome_second_order_sleep.parameter_connection_between_region,
                                  parameter_76connectome_second_order_sleep.parameter_coupling,
                                  parameter_76connectome_second_order_sleep.parameter_integrator,
                                  parameter_76connectome_second_order_sleep.parameter_monitor)
# tools_simulation.run_simulation(simulator,
#                                 10000.0,
#                                 parameter_76connectome_second_order_sleep.parameter_simulation,
#                                 parameter_76connectome_second_order_sleep.parameter_monitor)

# tools_simulation.print_all(parameter_76connectome_second_order_sleep.parameter_simulation['path_result'],
#                            0.0,10000.0,0,0)
#
# tools_simulation.print_all(parameter_76connectome_second_order_sleep.parameter_simulation['path_result'],
#                            0.0,10000.0,
#                            0,5)
#
# tools_simulation.print_EI_one(parameter_76connectome_second_order_sleep.parameter_simulation['path_result'],
#                            0.0,10000.0,
#                            0,0)

