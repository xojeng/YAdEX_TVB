import parameter_oneRegion_second_order_Wake
import tools_simulation
import numpy as np

parameter_oneRegion_second_order_Wake.parameter_simulation['path_result'] = './test/'



# simulator = tools_simulation.init(parameter_oneRegion_second_order_Wake.parameter_simulation,
#                                   parameter_oneRegion_second_order_Wake.parameter_model,
#                                   parameter_oneRegion_second_order_Wake.parameter_connection_between_region,
#                                   parameter_oneRegion_second_order_Wake.parameter_coupling,
#                                   parameter_oneRegion_second_order_Wake.parameter_integrator,
#                                   parameter_oneRegion_second_order_Wake.parameter_monitor)
#
# tools_simulation.run_simulation(simulator,
#                                 30000.0,
#                                 parameter_oneRegion_second_order_Wake.parameter_simulation,
#                                 parameter_oneRegion_second_order_Wake.parameter_monitor)
# #
# tools_simulation.print_EI_one(parameter_oneRegion_second_order_Wake.parameter_simulation['path_result'],
#                            0.0,30000.0,
#                            0,0)
#
# tools_simulation.print_bistability(parameter_oneRegion_second_order_Wake.parameter_model)

# tools_simulation.print_bistability(parameter_one_second_order.parameter_model)

result = tools_simulation.get_result(parameter_oneRegion_second_order_Wake.parameter_simulation['path_result'],0.0,10000.0)
time_raw = result[0][0]
raw = result[0][1]
time_bold = result[1][0]
bold_result = result[1][1]
import bold
monitor = bold.Bold()
monitor.period = parameter_oneRegion_second_order_Wake.parameter_integrator['dt']*2000.0
monitor.config_for_sim(parameter_oneRegion_second_order_Wake.parameter_connection_between_region['number_of_regions'],
                       parameter_oneRegion_second_order_Wake.parameter_integrator['dt']
                       )
save_time = []
save_bold =[]
for i in np.arange(1,result[0][0].shape[0],1):
    b = monitor.sample(i+1,raw[i][0])
    if b is not None:
        save_time.append(b[0])
        save_bold.append(b[1])

if any(save_time != time_bold):
    print('bad ')

if np.max(np.abs(save_bold - np.squeeze(bold_result,1))) > 5e-5:
    print('bad ')