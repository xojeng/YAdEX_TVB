parameter_simulation={
    'path_result':'./test_simulation',
    'seed':10, # the seed for the random generator
    'save_time': 10000.0, # the time of simulation in each file
}

parameter_model ={
    #order of the model
    'order':1,
    #parameter of the model
    'g_L':10.0,
    'E_L_e':-65.0,
    'E_L_i':-65.0,
    'C_m':200.0,
    'b_e':60.0,
    'a_e':4.0,
    'b_i':0.0,
    'a_i':0.0,
    'tau_w_e':500.0,
    'tau_w_i':1.0,
    'E_e':0.0,
    'E_i':0.0,
    'Q_e':1.5,
    'Q_i':5.0,
    'tau_e':5.0,
    'tau_i':5.0,
    'N_tot':10000,
    'p_connect':0.05,
    'g':0.2,
    'T':20.0,
    'P_e':[-0.04983106,  0.005063550882777035,  -0.023470121807314552, 0.0022951513725067503,
            -0.0004105302652029825,  0.010547051343547399,  -0.03659252821136933,
            0.007437487505797858,  0.001265064721846073, -0.04072161294490446],
    'P_i':[-0.05149122024209484,  0.004003689190271077, -0.008352013668528155,0.0002414237992765705,
                    -0.0005070645080016026,  0.0014345394104282397, -0.014686689498949967,
                    0.004502706285435741,0.0028472190352532454, -0.015357804594594548],
    'external_input_ex_ex':0.000,
    'external_input_ex_in':0.000,
    'external_input_in_ex':0.000,
    'external_input_in_in':0.000,
    #Initial condition :
    'initial_condition':{
        "E": [0.0,0.0],
        "I": [0.0,0.0],
        "W_e": [0.0,0.0],
        "W_i": [0.0,0.0]}
}

parameter_connection_between_region={
    ## CONNECTIVITY
    # connectivity by default
    'default':False,
    #from file (repertory with following files : tract_lengths.npy and weights.npy)
    'from_file':False,
    'path':'', #repertory of the files
    # File description
    'number_of_regions':2, # number of region
    # lenghts of tract between region : dimension => (number_of_regions, number_of_regions)
    'tract_lengths':[[1.0,1.0],[1.0,1.0]],
    # weight along the tract : dimension => (number_of_regions, number_of_regions)
    'weights':[[0.0,0.0],[0.0,0.0]],
    # speed of along long range connection
    'speed':4.0,
}

parameter_coupling={
    ##COUPLING
    'type':'Linear', # choice : Linear, Scaling, HyperbolicTangent, Sigmoidal, SigmoidalJansenRit, PreSigmoidal, Difference, Kuramoto
    'parameter':{'a':0.0,
                 'b':0.0}
}

parameter_integrator={
    ## INTEGRATOR
    'type':'Heun', # choice : Heun, Euler
    'stochastic':True,
    'noise_type': 'Ornstein_Ulhenbeck_process', # choice : Additive,Multiplicative,
    'noise_parameter':{
        'nsig':[5.0e-6,5.0e-6,0.0,0.0],
        'ntau':0.0,
        'dt':0.5
                        },
    'dt':0.5 # in ms
}

parameter_monitor= {
    'Raw':True,
    'TemporalAverage':True,
    'parameter_TemporalAverage':{
        'variables_of_interest':[0,1,2,3],
        'period':parameter_integrator['dt']*10.0
    },
    'Bold':False,
    'parameter_Bold':{
        'variables_of_interest':[0],
        'period':parameter_integrator['dt']*2000.0
    }
    # More monitor can be added
}