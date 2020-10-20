import parameter_default_connectome_second_order
import tools_simulation
import numpy as np

simulator = tools_simulation.init(parameter_default_connectome_second_order.parameter_simulation,
                                  parameter_default_connectome_second_order.parameter_model,
                                  parameter_default_connectome_second_order.parameter_connection_between_region,
                                  parameter_default_connectome_second_order.parameter_coupling,
                                  parameter_default_connectome_second_order.parameter_integrator,
                                  parameter_default_connectome_second_order.parameter_monitor)

tools_simulation.run_simulation(simulator,
                                1200000.0,
                                parameter_default_connectome_second_order.parameter_simulation,
                                parameter_default_connectome_second_order.parameter_monitor)

