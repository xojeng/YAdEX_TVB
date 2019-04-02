import parameter_one_first_order
import tools_simulation

simulator = tools_simulation.init(parameter_one_first_order.parameter_simulation,
                                  parameter_one_first_order.parameter_model,
                                  parameter_one_first_order.parameter_connection_between_region,
                                  parameter_one_first_order.parameter_coupling,
                                  parameter_one_first_order.parameter_integrator,
                                  parameter_one_first_order.parameter_monitor)

tools_simulation.run_simulation(simulator,
                                20000.0,
                                parameter_one_first_order.parameter_simulation,
                                parameter_one_first_order.parameter_monitor)