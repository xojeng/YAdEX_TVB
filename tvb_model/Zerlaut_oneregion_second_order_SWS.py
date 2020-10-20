import parameter_oneRegion_second_order_SlowWaves
import tools_simulation
import numpy as np

parameter_oneRegion_second_order_SlowWaves.parameter_simulation['path_result'] = './sleep_TA2/'
# parameter_one_second_order.parameter_model['b_e']=1.0
# parameter_one_second_order.parameter_model['a_e']=0.0
# parameter_one_second_order.parameter_model['tau_w_e']=500.0
# parameter_one_second_order.parameter_model['P_e']=[-0.0498,  0.00506,  -0.025, 0.0014, -0.00041,  0.0105,  -0.036, 0.0074,  0.0012, -0.0407]
# parameter_one_second_order.parameter_model['P_i']=[-0.0514,  0.004, -0.0083,0.0002,-0.0005,  0.0014, -0.0146, 0.0045,0.0028, -0.0153]
# parameter_one_second_order.parameter_model['Q_e']=1.0
# parameter_one_second_order.parameter_model['E_L_e']=-50.0
# parameter_one_second_order.parameter_model['E_L_i']=-70.0
# # parameter_one_second_order.parameter_model['a_e']=0.0
# # parameter_one_second_order.parameter_model['T']=40.0
# parameter_one_second_order.parameter_integrator['noise_parameter']['nsig']=[1.0e-8,1.0e-8,0.0,0.0,0.0,0.0,0.0,]
# parameter_one_second_order.parameter_model['initial_condition']={
#          "E": [0.01,0.01],
#          "I": [0.01,0.01],
#          "C_ee": [0.0,0.0],
#          "C_ei": [0.0,0.0],
#          "C_ii": [0.0,0.0],
#          "W_e": [0.0,0.0],
#          "W_i": [0.0,0.0]}



simulator = tools_simulation.init(parameter_oneRegion_second_order_SlowWaves.parameter_simulation,
                                  parameter_oneRegion_second_order_SlowWaves.parameter_model,
                                  parameter_oneRegion_second_order_SlowWaves.parameter_connection_between_region,
                                  parameter_oneRegion_second_order_SlowWaves.parameter_coupling,
                                  parameter_oneRegion_second_order_SlowWaves.parameter_integrator,
                                  parameter_oneRegion_second_order_SlowWaves.parameter_monitor)

tools_simulation.run_simulation(simulator,
                                30000.0,
                                parameter_oneRegion_second_order_SlowWaves.parameter_simulation,
                                parameter_oneRegion_second_order_SlowWaves.parameter_monitor)

tools_simulation.print_EI_one(parameter_oneRegion_second_order_SlowWaves.parameter_simulation['path_result'],
                           0.0,30000.0,
                           0,0)
#
tools_simulation.print_bistability(parameter_oneRegion_second_order_SlowWaves.parameter_model)

# tools_simulation.print_bistability(parameter_one_second_order.parameter_model)