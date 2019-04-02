import Zerlaut
import tvb.simulator.lab as lab
import numpy.random as rgn
import numpy as np
import noise as my_noise
import json
import os

def init(parameter_simulation,parameter_model,parameter_connection_between_region,parameter_coupling,
         parameter_integrator,parameter_monitor):
    '''
    Initialise the simulator with parameter

    :param parameter_simulation: parameters for the simulation
    :param parameter_model: parameters for the model
    :param parameter_connection_between_region: parameters for the connection between nodes
    :param parameter_coupling: parameters for the coupling of equations
    :param parameter_integrator: parameters for the intergator of the equation
    :param parameter_monitor: parameters for the monitors
    :return: the simulator initialize
    '''
    ## initialise the random generator
    rgn.seed(parameter_simulation['seed'])

    ## Model
    if parameter_model['order'] == 1:
        model = Zerlaut.Zerlaut_adaptation_first_order(variables_of_interest='E I W_e W_i'.split())
    elif parameter_model['order'] == 2:
        model = Zerlaut.Zerlaut_adaptation_second_order(variables_of_interest='E I C_ee C_ei C_ii W_e W_i'.split())
    else:
        raise Exception('Bad order for the model')

    model.g_L=parameter_model['g_L']
    model.E_L_e=parameter_model['E_L_e']
    model.E_L_i=parameter_model['E_L_i']
    model.C_m=parameter_model['C_m']
    model.b_e=parameter_model['b_e']
    model.a_e=parameter_model['a_e']
    model.b_i=parameter_model['b_i']
    model.a_i=parameter_model['a_i']
    model.tau_w_e=parameter_model['tau_w_e']
    model.tau_w_i=parameter_model['tau_w_i']
    model.E_e=parameter_model['E_e']
    model.E_i=parameter_model['E_i']
    model.Q_e=parameter_model['Q_e']
    model.Q_i=parameter_model['Q_i']
    model.tau_e=parameter_model['tau_e']
    model.tau_i=parameter_model['tau_i']
    model.N_tot=parameter_model['N_tot']
    model.p_connect=parameter_model['p_connect']
    model.g=parameter_model['g']
    model.T=parameter_model['T']
    model.P_e=parameter_model['P_e']
    model.P_i=parameter_model['P_i']
    model.external_input_ex_ex=parameter_model['external_input_ex_ex']
    model.external_input_ex_in=parameter_model['external_input_ex_in']
    model.external_input_in_ex=parameter_model['external_input_in_ex']
    model.external_input_in_in=parameter_model['external_input_in_in']
    model.state_variable_range['E'] = parameter_model['initial_condition']['E']
    model.state_variable_range['I'] = parameter_model['initial_condition']['I']
    if parameter_model['order'] == 2:
        model.state_variable_range['C_ee'] = parameter_model['initial_condition']['C_ee']
        model.state_variable_range['C_ei'] = parameter_model['initial_condition']['C_ei']
        model.state_variable_range['C_ii'] = parameter_model['initial_condition']['C_ii']
    model.state_variable_range['W_e'] = parameter_model['initial_condition']['W_e']
    model.state_variable_range['W_i'] = parameter_model['initial_condition']['W_i']

    ## Connection
    if parameter_connection_between_region['default']:
        connection = lab.connectivity.Connectivity(load_default=True)
    elif parameter_connection_between_region['from_file']:
        path = parameter_connection_between_region['path']
        tract_lengths = np.load(path+'/tract_lengths.npy')
        weights = np.load(path+'/weights.npy')
        connection = lab.connectivity.Connectivity(number_of_regions=tract_lengths.shape[0],
                                               tract_lengths=tract_lengths,
                                               weights=weights)
    else:
        connection = lab.connectivity.Connectivity(number_of_regions=parameter_connection_between_region['number_of_regions'],
                                               tract_lengths=parameter_connection_between_region['tract_lengths'],
                                               weights=parameter_connection_between_region['weights'],)
    connection.speed = parameter_connection_between_region['speed']

    ## Coupling
    if parameter_coupling['type'] == 'Linear':
        coupling = lab.coupling.Linear(a=parameter_coupling['parameter']['a'],
                                       b=parameter_coupling['parameter']['b'])
    elif parameter_coupling['type'] == 'Scaling':
        coupling = lab.coupling.Scaling(a=parameter_coupling['parameter']['a'])
    elif parameter_coupling['type'] == 'HyperbolicTangent':
        coupling = lab.coupling.HyperbolicTangent(a=parameter_coupling['parameter']['a'],
                                       b=parameter_coupling['parameter']['b'],
                                       midpoint=parameter_coupling['parameter']['midpoint'],
                                       sigma= parameter_coupling['parameter']['sigma'],)
    elif parameter_coupling['type'] == 'Sigmoidal':
        coupling = lab.coupling.Sigmoidal(a=parameter_coupling['parameter']['a'],                                       b=parameter_coupling['b'],
                                       midpoint=parameter_coupling['parameter']['midpoint'],
                                       sigma= parameter_coupling['parameter']['sigma'],
                                       cmin=parameter_coupling['parameter']['cmin'],
                                       cmax=parameter_coupling['parameter']['cmax'])
    elif parameter_coupling['type'] == 'SigmoidalJansenRit':
        coupling = lab.coupling.SigmoidalJansenRit(a=parameter_coupling['parameter']['a'],                                       b=parameter_coupling['b'],
                                       midpoint=parameter_coupling['parameter']['midpoint'],
                                       r= parameter_coupling['parameter']['r'],
                                       cmin=parameter_coupling['parameter']['cmin'],
                                       cmax=parameter_coupling['parameter']['cmax'])
    elif parameter_coupling['type'] == 'PreSigmoidal':
        coupling = lab.coupling.PreSigmoidal(H=parameter_coupling['parameter']['H'],                                       b=parameter_coupling['b'],
                                       Q=parameter_coupling['parameter']['Q'],
                                       G= parameter_coupling['parameter']['G'],
                                       P=parameter_coupling['parameter']['P'],
                                       theta=parameter_coupling['parameter']['theta'],
                                       dynamic=parameter_coupling['parameter']['dynamic'],
                                       globalT=parameter_coupling['parameter']['globalT'],
                                             )
    elif parameter_coupling['type'] == 'Difference':
        coupling = lab.coupling.Difference(a=parameter_coupling['parameter']['a'])
    elif parameter_coupling['type'] == 'Kuramoto':
        coupling = lab.coupling.Kuramoto(a=parameter_coupling['parameter']['a'])
    else:
        raise Exception('Bad type for the coupling')

    ## Integrator
    if not parameter_integrator['stochastic']:
        if parameter_integrator['type'] == 'Heun':
            integrator = lab.integrators.HeunDeterministic(dt=parameter_integrator['dt'])
        elif parameter_integrator['type'] == 'Euler':
             integrator = lab.integrators.EulerDeterministic(dt=parameter_integrator['dt'])
        else:
            raise Exception('Bad type for the integrator')
    else:
        if parameter_integrator['noise_type'] == 'Ornstein_Ulhenbeck_process':
            noise = my_noise.Ornstein_Ulhenbeck_process(nsig=parameter_integrator['noise_parameter']['nsig'],
                                                        ntau=parameter_integrator['noise_parameter']['ntau'],)
        elif parameter_integrator['noise_type'] == 'Additive':
            noise = lab.noise.Additive(nsig=parameter_integrator['noise_parameter']['nsig'],
                                                        ntau=parameter_integrator['noise_parameter']['ntau'],)
        elif parameter_integrator['noise_type'] == 'Ornstein_Ulhenbeck_process':
            noise = lab.noise.Multiplicative(nsig=parameter_integrator['noise_parameter']['nsig'],
                                                        ntau=parameter_integrator['noise_parameter']['ntau'],
                                                        b=parameter_integrator['noise_parameter']['b'],
                                                        )
        else:
            raise Exception('Bad type for the noise')

        if parameter_integrator['type'] == 'Heun':
            integrator = lab.integrators.HeunStochastic(noise=noise,dt=parameter_integrator['dt'])
        elif parameter_integrator['type'] == 'Euler':
             integrator = lab.integrators.EulerStochastic(noise=noise,dt=parameter_integrator['dt'])
        else:
            raise Exception('Bad type for the integrator')

    ## Monitors
    monitors =[]
    if parameter_monitor['Raw']:
        monitors.append(lab.monitors.Raw())
    if parameter_monitor['TemporalAverage']:
        monitor_TAVG = lab.monitors.TemporalAverage(
            variables_of_interest=parameter_monitor['parameter_TemporalAverage']['variables_of_interest'],
            period=parameter_monitor['parameter_TemporalAverage']['period'])
        monitors.append(monitor_TAVG)
    if parameter_monitor['Bold']:
        monitor_Bold = lab.monitors.Bold(
            variables_of_interest=parameter_monitor['parameter_Bold']['variables_of_interest'],
            period=parameter_monitor['parameter_Bold']['period'])
        monitors.append(monitor_Bold)

    #save the parameters in on file
    if not os.path.exists(parameter_simulation['path_result']):
        os.mkdir(parameter_simulation['path_result'])
    f = open(parameter_simulation['path_result']+'/parameter.py',"w")
    for name,dic in [('parameter_simulation',parameter_simulation),
                     ('parameter_model',parameter_model),
                     ('parameter_connection_between_region',parameter_connection_between_region),
                     ('parameter_coupling',parameter_coupling),
                     ('parameter_integrator',parameter_integrator),
                     ('parameter_monitor',parameter_monitor)]:
        f.write(name+' = ')
        json.dump(dic, f)
        f.write("\n")
    f.close()

    #initialize the simulator
    simulator = lab.simulator.Simulator(model = model, connectivity = connection,
                          coupling = coupling, integrator = integrator, monitors = monitors)
    simulator.configure()
    return simulator

def run_simulation(simulator, time, parameter_simulation,parameter_monitor):
    '''
    run a simulation
    :param simulator: the simulator already initialize
    :param time: the time of simulation
    :param parameter_simulation: the parameter for the simulation
    :param parameter_monitor: the parameter for the monitor
    '''
    # save the initial condition
    np.save(parameter_simulation['path_result']+'/step_init.npy',simulator.history.buffer)
    # check how many monitor it's used
    nb_monitor = parameter_monitor['Raw'] + parameter_monitor['TemporalAverage'] + parameter_monitor['Bold']
    # initialise the variable for the saving the result
    save_result =[]
    for i in range(nb_monitor):
        save_result.append([])
    # run the simulation
    count = 0
    for result in simulator(simulation_length=time):
        for i in range(nb_monitor):
            if result[i] is not None:
                save_result[i].append(result[i])
        #save the result in file
        if result[0][0] >= parameter_simulation['save_time']*(count+1): #check if the time for saving at some time step
            print('simulation time :'+str(result[0][0])+'\r')
            np.save(parameter_simulation['path_result']+'/step_'+str(count)+'.npy',save_result)
            save_result =[]
            for i in range(nb_monitor):
                save_result.append([])
            count +=1
    # save the last part
    np.save(parameter_simulation['path_result']+'/step_'+str(count)+'.npy',save_result)

