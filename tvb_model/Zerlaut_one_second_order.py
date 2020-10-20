import parameter_default_connectome_second_order as parameter
import tools_simulation
import numpy as np
import time

simulator = tools_simulation.init(parameter.parameter_simulation,
                                  parameter.parameter_model,
                                  parameter.parameter_connection_between_region,
                                  parameter.parameter_coupling,
                                  parameter.parameter_integrator,
                                  parameter.parameter_monitor)
print(time.ctime())
tools_simulation.run_simulation(simulator,
                                10000.0,
                                parameter.parameter_simulation,
                                parameter.parameter_monitor)

print(time.ctime())
# tools_simulation.get_result(parameter.parameter_simulation['path_result'],0.0,1000.0)
