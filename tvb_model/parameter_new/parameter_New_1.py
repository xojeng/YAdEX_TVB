parameter_simulation={
    'path_result':'./simulation_default',
    'seed':10, # the seed for the random generator
    'save_time': 10000.0, # the time of simulation in each file
}

parameter_model ={
    #order of the model
    'order':2,
    #parameter of the model
    'g_L':10.0,
    'E_L_e':-67.0,
    'E_L_i':-65.0,
    'C_m':200.0,
    'b_e':10.0,
    'a_e':0.0,
    'b_i':0.0,
    'a_i':0.0,
    'tau_w_e':500.0,
    'tau_w_i':1.0,
    'E_e':0.0,
    'E_i':-80.0,
    'Q_e':1.0,
    'Q_i':6.5,
    'tau_e':5.0,
    'tau_i':5.0,
    'N_tot':10000,
    'p_connect_e':0.10,
    'p_connect_i':0.095,# True value : 0.05,
    'g':0.2,
    'T':20.0,
    'P_e':[-0.05514735,  0.00444736,  0.00079013,  0.01183319,  0.00145571,  0.00155584, -0.00679113,  0.00091574, -0.00475571,  0.00332159],
    'P_i':[-0.05637474,  0.005188,   -0.00342982,  0.01282597, -0.00142761, -0.00201966, -0.02266407,  0.00330289,  0.00486843, -0.00908176],
    'external_input_ex_ex':0.000,
    'external_input_ex_in':0.000,
    'external_input_in_ex':0.000,
    'external_input_in_in':0.000,
    'K_ext_e':400,
    'K_ext_i':0,
    #Initial condition :
    'initial_condition':{
        "E": [0.001, 0.0001],"I": [0.001, 0.0001],"C_ee": [0.0,0.0],"C_ei": [0.0,0.0],"C_ii": [0.0,0.0],"W_e": [0.0, 0.0],"W_i": [0.0,0.0]}
}

parameter_connection_between_region={
    ## CONNECTIVITY
    # connectivity by default
    'default':False,
    #from file (repertory with following files : tract_lengths.npy and weights.npy)
    'from_file':False,
    'from_h5':True,
    'path':'/home/kusch/Documents/project/Zerlaut/YAdEX_TVB/data/hcp-001.zip', #repertory of the files
    # File description
    'number_of_regions':0, # number of region
    # lenghts of tract between region : dimension => (number_of_regions, number_of_regions)
    'tract_lengths':[],
    # weight along the tract : dimension => (number_of_regions, number_of_regions)
    'weights':[],
    # speed of along long range connection
    'speed':4.0,
    'normalised':True
}

parameter_coupling={
    ##COUPLING
    'type':'Linear', # choice : Linear, Scaling, HyperbolicTangent, Sigmoidal, SigmoidalJansenRit, PreSigmoidal, Difference, Kuramoto
    'parameter':{'a':5.0e-1,
                 'b':0.0}
}

# parameter_integrator={
#     ## INTEGRATOR
#     'type':'RungeKutta4thOrderDeterministic', # choice : Heun, Euler
#     'dt': 0.5
# }

parameter_integrator={
    ## INTEGRATOR
    'type':'Heun', # choice : Heun, Euler
    'stochastic':True,
    'noise_type': 'Ornstein_Ulhenbeck_process', # choice : Additive,Multiplicative,
    'noise_parameter':{
        'mu':[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], # in KHz, TAedit: decreased value from 1.0e-10
        'nsig':[0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], # in KHz, TAedit: decreased value from 1.0e-10
        'tau_OU':0.1,
        'ntau':0.0,
        'weight':[1.e-2]
                        },
    'dt': 0.5 # in ms
}

parameter_monitor= {
    'Raw':True,
    'TemporalAverage':True,
    'parameter_TemporalAverage':{
        'variables_of_interest':[0,1,2,3,4,5,6],
        'period':parameter_integrator['dt']*10.0
    },
    'Bold':False,
    'parameter_Bold':{
        'variables_of_interest':[0],
        'period':parameter_integrator['dt']*2000.0
    }
    # More monitor can be added
}