import tools_simulation
import parameter_Matteo
import numpy as np
import utils
import matplotlib.pyplot as plt

t=20000.0

bvals = np.linspace(0,40,10)

for simnum in range(len(bvals)):

    invals = np.linspace(0.00005, 0.00015, 10)

    for input in range(len(invals)):

        parameter_Matteo.parameter_model['b_e'] = bvals[simnum]
        print(parameter_Matteo.parameter_coupling.keys())
        parameter_Matteo.parameter_coupling['parameter']['a'] = 1.0
        #parameter_Matteo.parameter_integrator['noise_parameter']['nsig'] = [invals[input], invals[input], 0.0, 0.0, 0.0, 0.0, 0.0]
        parameter_Matteo.parameter_model['external_input_in_ex'] = invals[input]
        parameter_Matteo.parameter_model['external_input_in_in'] = invals[input]
        parameter_Matteo.parameter_simulation['path_result']='/mnt/usb-WD_easystore_25FC_575856314136373548355350-0:0-part1/Less_noise_v07/b_'+str(bvals[simnum])+'input_'+str(invals[input])
        print(parameter_Matteo.parameter_simulation['path_result'])

        Nregions = 84
        parameter_stimulation = {}
        parameter_stimulation['onset'] = 100000 # onset time in ms    `1
        parameter_stimulation['tau'] = 10 # stimulus duration in ms
        parameter_stimulation['T'] = 100000 # interstimulus interval in ms
        parameter_stimulation['weights'] = np.zeros(Nregions)
        parameter_stimulation['weights'][59] = 0 #Previous : 1e-3

        simulator = tools_simulation.init(parameter_Matteo.parameter_simulation,
                                         parameter_Matteo.parameter_model,
                                         parameter_Matteo.parameter_connection_between_region,
                                         parameter_Matteo.parameter_coupling,
                                         parameter_Matteo.parameter_integrator,
                                         parameter_Matteo.parameter_monitor,
                                         parameter_stimulation = parameter_stimulation)

        tools_simulation.run_simulation(simulator,
                                        t,
                                        parameter_Matteo.parameter_simulation,
                                        parameter_Matteo.parameter_monitor)

# tools_simulation.print_all(parameter_Matteo.parameter_simulation['path_result'],
#                            0.0,t,
#                            0,simulator.connectivity)

# # tools_simulation.print_all(parameter_Matteo.parameter_simulation['path_result'],
# #                            1.0,t,
# #                            0,simulator.connectivity)
#
# tools_simulation.print_EI_one(parameter_Matteo.parameter_simulation['path_result'],
#                              0.0,t,
#                              0,59)
# # #
# # the_data = tools_simulation.get_result(parameter_Matteo.parameter_simulation['path_result'], 0.0, t)
# # np.mean(the_data,0)

# # plot on brains keep this forever
# #plot change in firing rate on brains
# the_data_pre = tools_simulation.get_result(parameter_Matteo.parameter_simulation['path_result'], 8700.0, 9000)[0][1]
# mean_pre = np.mean(the_data_pre,0)[0,:]
# the_data_post = tools_simulation.get_result(parameter_Matteo.parameter_simulation['path_result'], 9001.0, 9301.0)[0][1]
# mean_post = np.mean(the_data_post,0)[0,:]
# the_data = mean_post/mean_pre
# remapping = np.array([np.nan,  8.,  9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19.,
#        20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 31., 32.,
#        33., 34., 35., 36., 37., 38., 39., 40., 41., np.nan, 50., 51., 52.,
#        53., 54., 55., 56., 57., 58., 59., 60., 61., 62., 63., 64., 65.,
#        66., 67., 68., 69., 70., 71., 72., 73., 74., 75., 76., 77., 78.,
#        79., 80., 81., 82., 83., np.nan,  0.,  1.,  2.,  3.,  4.,  5.,  6.,
#         7., 42., 43., 44., 45., 46., 47., 48., 49.])
#
# other=np.array([71, 72, 73, 74, 75, 76, 77, 78,  1,  2,  3,  4,  5,  6,  7,  8,  9,
#        10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
#        27, 28, 29, 30, 31, 32, 33, 34, 79, 80, 81, 82, 83, 84, 85, 86, 36,
#        37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,
#        54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69])
# new_data = np.zeros((remapping.shape[0]+1,))
# new_data[other]  = the_data
# utils.multiview(new_data, shaded=False)
# plt.show()
# print(the_data)

# for time_step in range(int(t*10)):
#     data  = the_data[0][1][time_step,0,utils.cortex.region_mapping]
#     utils.multiview(data, zmin=0.0, zmax = 0.02, shaded=False)
#     # plt.show()
#     plt.savefig(parameter_Matteo.parameter_simulation['path_result']+'/time_step_'+str(time_step)+'.png')
#     plt.close('all')

