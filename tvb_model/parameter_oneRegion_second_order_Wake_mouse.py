parameter_simulation={
    'path_result':'/home/mansi/summer_project_2019/zerlaut/data/mouse/',
    'seed':10, # the seed for the random generator
    'save_time': 2000.0, # the time of simulation in each file
}

parameter_model ={
    #order of the model
    'order':2,
    #parameter of the model
    'g_L':10.0,
    'E_L_e':-50.0,
    'E_L_i':-70.0,
    'C_m':200.0,
    'b_e':1.0,
    'a_e':0.0,
    'b_i':0.0,
    'a_i':0.0,
    'tau_w_e':500.0,
    'tau_w_i':1.0,
    'E_e':0.0,
    'E_i':-80.0,
    'Q_e':1.0,
    'Q_i':5.0,
    'tau_e':5.0,
    'tau_i':5.0,
    'N_tot':10000,
    'p_connect':0.05,
    'g':0.2,
    'T':20.0,
    'P_e':[-0.0498, 0.00506, -0.025, 0.0014, -0.00041, 0.0105, -0.036, 0.0074, 0.0012, -0.0407],
    'P_i':[-0.0514, 0.004, -0.0083, 0.0002, -0.0005, 0.0014, -0.0146, 0.0045, 0.0028, -0.0153],
    'external_input_ex_ex':0.000,
    'external_input_ex_in':0.000,
    'external_input_in_ex':0.000,
    'external_input_in_in':0.000,
    #Initial condition :
    'initial_condition':{"E": [0.01, 0.01], "I": [0.01, 0.01], "C_ii": [0.0, 0.0], "W_e": [0.0, 0.0], "C_ee": [0.0, 0.0], "C_ei": [0.0, 0.0], "W_i": [0.0, 0.0]}
}

parameter_connection_between_region={
    ## CONNECTIVITY
    # connectivity by default
    'default':False,
    #from file (repertory with following files : tract_lengths.npy and weights.npy)
    'from_file':False,
    'from_h5':True,
    'path':'/home/mansi/summer_project_2019/TVB_Distribution/tvb_data/lib/python2.7/site-packages/tvb_data/mouse/allen_2mm/Connectivity.h5', #repertory of the files
    # File description
    'number_of_regions':98, # number of region
    # lenghts of tract between region : dimension => (number_of_regions, number_of_regions)
    'tract_lengths':[],
    # weight along the tract : dimension => (number_of_regions, number_of_regions)
    'weights':[],
    # speed of along long range connection
    'speed':3.0,
}

parameter_coupling={
    ##COUPLING
    'type':'Linear', # choice : Linear, Scaling, HyperbolicTangent, Sigmoidal, SigmoidalJansenRit, PreSigmoidal, Difference, Kuramoto
    'parameter':{'a':0.096,
                 'b':0.0}
}

parameter_integrator={
    ## INTEGRATOR
    'type':'Heun', # choice : Heun, Euler
    'stochastic':True,
    'noise_type': 'Ornstein_Ulhenbeck_process', # choice : Additive,Multiplicative,
    'noise_parameter':{
        'nsig':[1e-08, 1e-08, 0.0, 0.0, 0.0, 0.0, 0.0],
        'ntau':0.0,
        'dt':0.1
                        },
    'dt':0.1 # in ms
}

parameter_monitor= {
    'Raw':True,
    'TemporalAverage':False,
    'parameter_TemporalAverage':{
        'variables_of_interest':[0,1,2,3],
        'period':parameter_integrator['dt']*10.0
    },
    'Bold':False,
    'parameter_Bold':{
        'variables_of_interest':[0],
        'period':parameter_integrator['dt']*20000.0
    }
    # More monitor can be added
}
