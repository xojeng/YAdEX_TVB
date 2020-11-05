from tvb_model_reference.simulation_file.compare_result.parameter_Matteo_original import Parameter as Parameter_second_order
from tvb_model_reference.simulation_file.compare_result.parameter_Matteo_original_first_order import Parameter as Parameter_first_order
import tvb_model_reference.src.tools_simulation as tools
import matplotlib.pyplot as plt
import copy


M_original = Parameter_first_order()
M_original.parameter_simulation['path_result']="./result/M_original"
M_mono = copy.deepcopy(M_original)
M_mono.parameter_model['E_L_e']=-67
M_mono.parameter_simulation['path_result']="./result/M_mono"
M_bi = copy.deepcopy(M_original)
M_bi.parameter_model['E_L_e']=-63
M_bi.parameter_simulation['path_result']="./result/M_bi"
TVB_original = copy.deepcopy(M_original)
TVB_original.parameter_model['matteo']=False
TVB_original.parameter_simulation['path_result']="./result/TVB_original"
TVB_mono = copy.deepcopy(TVB_original)
TVB_mono.parameter_model['E_L_e']=-67
TVB_mono.parameter_simulation['path_result']="./result/TVB_mono"
TVB_bi = copy.deepcopy(TVB_original)
TVB_bi.parameter_model['E_L_e']=-63
TVB_bi.parameter_simulation['path_result']="./result/TVB_bi"
run_sim = 500.0
for parameters in [M_original,M_bi,M_mono,TVB_original,TVB_bi,TVB_mono]:
    plt.figure()
    simulator = tools.init(parameters.parameter_simulation,
                                  parameters.parameter_model,
                                  parameters.parameter_connection_between_region,
                                  parameters.parameter_coupling,
                                  parameters.parameter_integrator,
                                  parameters.parameter_monitor)
    tools.run_simulation(simulator,
                                run_sim,
                                parameters.parameter_simulation,
                                parameters.parameter_monitor)
    plt.title("second order "+ parameters.parameter_simulation['path_result'])
    tools.print_EI_one(parameters.parameter_simulation['path_result'],
                           0.0,run_sim,0,0)

M_original = Parameter_second_order()
M_original.parameter_simulation['path_result']="./result/M_original"
M_mono = copy.deepcopy(M_original)
M_mono.parameter_model['E_L_e']=-67
M_mono.parameter_simulation['path_result']="./result/M_mono"
M_bi = copy.deepcopy(M_original)
M_bi.parameter_model['E_L_e']=-63
M_bi.parameter_simulation['path_result']="./result/M_bi"
TVB_original = copy.deepcopy(M_original)
TVB_original.parameter_model['matteo']=False
TVB_original.parameter_simulation['path_result']="./result/TVB_original"
TVB_mono = copy.deepcopy(TVB_original)
TVB_mono.parameter_model['E_L_e']=-67
TVB_mono.parameter_simulation['path_result']="./result/TVB_mono"
TVB_bi = copy.deepcopy(TVB_original)
TVB_bi.parameter_model['E_L_e']=-63
TVB_bi.parameter_simulation['path_result']="./result/TVB_bi"
run_sim = 500.0
for parameters in [M_original,M_bi,M_mono,TVB_original,TVB_bi,TVB_mono]:
    simulator = tools.init(parameters.parameter_simulation,
                                  parameters.parameter_model,
                                  parameters.parameter_connection_between_region,
                                  parameters.parameter_coupling,
                                  parameters.parameter_integrator,
                                  parameters.parameter_monitor)
    tools.run_simulation(simulator,
                                run_sim,
                                parameters.parameter_simulation,
                                parameters.parameter_monitor)
    plt.title("first order "+ parameters.parameter_simulation['path_result'])
    tools.print_EI_one(parameters.parameter_simulation['path_result'],
                           0.0,run_sim,0,0)
plt.show()