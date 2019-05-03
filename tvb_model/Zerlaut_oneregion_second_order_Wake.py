import parameter_oneRegion_second_order_Wake
import tools_simulation
import numpy as np

parameter_oneRegion_second_order_Wake.parameter_simulation['path_result'] = './wake_TA2/'



simulator = tools_simulation.init(parameter_oneRegion_second_order_Wake.parameter_simulation,
                                  parameter_oneRegion_second_order_Wake.parameter_model,
                                  parameter_oneRegion_second_order_Wake.parameter_connection_between_region,
                                  parameter_oneRegion_second_order_Wake.parameter_coupling,
                                  parameter_oneRegion_second_order_Wake.parameter_integrator,
                                  parameter_oneRegion_second_order_Wake.parameter_monitor)

tools_simulation.run_simulation(simulator,
                                30000.0,
                                parameter_oneRegion_second_order_Wake.parameter_simulation,
                                parameter_oneRegion_second_order_Wake.parameter_monitor)

tools_simulation.print_EI_one(parameter_oneRegion_second_order_Wake.parameter_simulation['path_result'],
                           0.0,30000.0,
                           0,0)
#
tools_simulation.print_bistability(parameter_oneRegion_second_order_Wake.parameter_model)

# tools_simulation.print_bistability(parameter_one_second_order.parameter_model)