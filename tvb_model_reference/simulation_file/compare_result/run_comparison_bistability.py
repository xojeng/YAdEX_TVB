from tvb_model_reference.simulation_file.compare_result.parameter_Matteo_original import Parameter
import tvb_model_reference.src.tools_simulation as tools
import matplotlib.pyplot as plt
import copy

M_original = Parameter()
M_mono = copy.deepcopy(M_original)
M_mono.parameter_model['E_L_e']=-67
M_bi = copy.deepcopy(M_original)
M_bi.parameter_model['E_L_e']=-63
TVB_original = copy.deepcopy(M_original)
TVB_original.parameter_model['matteo']=False
TVB_mono = copy.deepcopy(TVB_original)
TVB_mono.parameter_model['E_L_e']=-67
TVB_bi = copy.deepcopy(TVB_original)
TVB_bi.parameter_model['E_L_e']=-63

data_M_original = tools.print_bistability(M_original.parameter_model,show=False)
data_M_mono = tools.print_bistability(M_mono.parameter_model,show=False)
data_M_bi = tools.print_bistability(M_bi.parameter_model,show=False)
data_TVB_original = tools.print_bistability(TVB_original.parameter_model,show=False)
data_TVB_mono = tools.print_bistability(TVB_mono.parameter_model,show=False)
data_TVB_bi = tools.print_bistability(TVB_bi.parameter_model,show=False)

plt.figure()
for name,data in [('matteo original',data_M_original),('matteo mono', data_M_mono), ('matteo bi',data_M_bi),('TVB original',data_TVB_original),('TVB mono',data_TVB_mono),('TVB bi',data_TVB_bi)]:
    plt.plot(data[0]*1e3,data[1]*1e3,label=name)
plt.plot(data[0]*1e3,data[0]*1e3,'k--')
plt.legend()
plt.show()
