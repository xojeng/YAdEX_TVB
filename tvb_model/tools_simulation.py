import Zerlaut
import tvb.simulator.lab as lab
import numpy.random as rgn
import numpy as np
import noise as my_noise
import json
import os
import subprocess
import sys
from scipy.optimize import fsolve

def init(parameter_simulation,parameter_model,parameter_connection_between_region,parameter_coupling,
         parameter_integrator,parameter_monitor, initial_condition=None):
    '''
    Initialise the simulator with parameter

    :param parameter_simulation: parameters for the simulation
    :param parameter_model: parameters for the model
    :param parameter_connection_between_region: parameters for the connection between nodes
    :param parameter_coupling: parameters for the coupling of equations
    :param parameter_integrator: parameters for the intergator of the equation
    :param parameter_monitor: parameters for the monitors
    :param initial_condition: the possibility to add an initial condition
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
        noise.random_stream.seed(parameter_simulation['seed'])

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
    subprocess.call(['./correct_parameter.sh',parameter_simulation['path_result']+'/parameter.py']) ##Warning can be do'nt find the script
    #initialize the simulator
    if initial_condition == None:
        simulator = lab.simulator.Simulator(model = model, connectivity = connection,
                          coupling = coupling, integrator = integrator, monitors = monitors)
    else:
        simulator = lab.simulator.Simulator(model = model, connectivity = connection,
                                            coupling = coupling, integrator = integrator,
                                            monitors = monitors, initial_conditions=initial_condition)
    simulator.configure()
    if initial_condition == None:
        # save the initial condition
        np.save(parameter_simulation['path_result']+'/step_init.npy',simulator.history.buffer)
    return simulator

def run_simulation(simulator, time, parameter_simulation,parameter_monitor):
    '''
    run a simulation
    :param simulator: the simulator already initialize
    :param time: the time of simulation
    :param parameter_simulation: the parameter for the simulation
    :param parameter_monitor: the parameter for the monitor
    '''
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

def get_result(path,time_begin,time_end):
    '''
    return the result of the simulation between the wanted time
    :param path: the folder of the simulation
    :param time_begin: the start time for the result
    :param time_end:  the ending time for the result
    :return: result of all monitor
    '''
    sys.path.append(path)
    from parameter import parameter_simulation,parameter_monitor
    sys.path.remove(path)
    count_begin = int(time_begin/parameter_simulation['save_time'])
    count_end = int(time_end/parameter_simulation['save_time'])+1
    nb_monitor = parameter_monitor['Raw'] + parameter_monitor['TemporalAverage'] + parameter_monitor['Bold']
    output =[]

    for count in range(count_begin,count_end):
        result = np.load(path+'/step_'+str(count)+'.npy',allow_pickle=True)
        for i in range(result.shape[0]):
            tmp = np.array(result[i])
            if len(tmp) != 0:
                tmp = tmp[np.where((time_begin <= tmp[:,0]) &  (tmp[:,0]<= time_end)),:]
                tmp_time = tmp[0][:,0]
                if tmp_time.shape[0] != 0:
                    one = tmp[0][:,1][0]
                    tmp_value = np.concatenate(tmp[0][:,1]).reshape(tmp_time.shape[0],one.shape[0],one.shape[1])
                    if len(output) == nb_monitor:
                        output[i]=[np.concatenate([output[i][0],tmp_time]),np.concatenate([output[i][1],tmp_value])]
                    else:
                        output.append([tmp_time,tmp_value])
    return output

def get_region(path,time_begin,time_end,region_nb):
    '''
    return the result of the simulation between the wanted time
    :param path: the folder of the simulation
    :param time_begin: the start time for the result
    :param time_end: the ending time for the result
    :param region_nb: interger of array of interger for the select the region of different regions
    :return: result of all monitor
    '''
    sys.path.append(path)
    from parameter import parameter_simulation,parameter_monitor
    sys.path.remove(path)
    count_begin = int(time_begin/parameter_simulation['save_time'])
    count_end = int(time_end/parameter_simulation['save_time'])+1
    nb_monitor = parameter_monitor['Raw'] + parameter_monitor['TemporalAverage'] + parameter_monitor['Bold']
    output =[]

    for count in range(count_begin,count_end):
        result = np.load(path+'/step_'+str(count)+'.npy')
        for i in range(result.shape[0]):
            tmp = np.array(result[i])
            tmp = tmp[np.where((time_begin <= tmp[:,0]) &  (tmp[:,0]<= time_end)),:]
            tmp_time = tmp[0][:,0]
            if tmp_time.shape[0] != 0:
                one = tmp[0][:,1][0]
                tmp_value = np.concatenate(tmp[0][:,1]).reshape(tmp_time.shape[0],one.shape[0],one.shape[1])
                tmp_value = tmp_value[:,:,region_nb]
                if len(output) == nb_monitor:
                    output[i]=[np.concatenate([output[i][0],tmp_time]),np.concatenate([output[i][1],tmp_value])]
                else:
                    output.append([tmp_time,tmp_value])
    return output

def print_all(path,time_begin,time_end,position_monitor,position_variable):
    '''
    print one value of on monitor for some times
    :param path: the folder of the simulation
    :param time_begin: the start time for the result
    :param time_end: the ending time for the result
    :param position_monitor: select the monitor
    :param position_variable: select the variable of monitor
    :return: nothing
    '''
    output = get_result(path,time_begin,time_end)
    import matplotlib.pyplot as plt
    plt.plot(output[position_monitor][0]*1e-3, # time in second
             output[position_monitor][1][:,position_variable,:]
             )
    plt.show()

def print_EI_one(path,time_begin,time_end,position_monitor,position_node):
    '''

    :param path: the folder of the simulation
    :param time_begin: the start time for the result
    :param time_end: the ending time for the result
    :param position_monitor: select the monitor
    :param position_node:
    :return: nothing
    '''
    output = get_result(path,time_begin,time_end)
    import matplotlib.pyplot as plt
    plt.plot(output[position_monitor][0]*1e-3, # time in second
             output[position_monitor][1][:,0,position_node],
             color='c')
    plt.plot(output[position_monitor][0]*1e-3, # time in second
             output[position_monitor][1][:,1,position_node],
             color='r')
    plt.figure()
    plt.plot(output[position_monitor][0]*1e-3, # time in second
             output[position_monitor][1][:,5,position_node],
             color='g')
    plt.show()

def print_region(path,time_begin,time_end,position_monitor,position_variable,nb_region):
    '''
    print one value of on monitor for some times
    :param path: the folder of the simulation
    :param time_begin: the start time for the result
    :param time_end: the ending time for the result
    :param position_monitor: select the monitor
    :param position_variable: select the variable of monitor
    :param nb_region: interger of array of interger for the select the region of different regions
    :return: nothing
    '''
    output = get_region(path,time_begin,time_end,nb_region)
    import matplotlib.pyplot as plt
    plt.plot(output[position_monitor][0]*1e-3, # time in second
             output[position_monitor][1][:,position_variable,:]
             )
    plt.show()

def print_bistability(parameter_model):
    '''
    print if the model is bistable or not
    :param parameter_model: parameters for the model
        (the parameter external_input_in_in and external_input_in_ex is taking in count)
    :return: nothing
    '''
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

    ##Solution equation
    def equation(x):
        return (model.TF_inhibitory((fezero+model.external_input_ex_in), x,0.0)-x)
        #return (model.TF_inhibitory((fezero+extinp), x,fezero*model.b_e*model.tau_w_e)-x) # with fix adaptation

    # raneg fo frequency in KHz
    xrange1 = np.arange(0.0, 200.0, 0.01)*1e-3
    feprimovec=0*xrange1
    fiprimovec=0*xrange1
    for i in range(len(xrange1)):
        fezero=xrange1[i]
        fizero=fsolve(equation, (fezero*2),xtol=1.e-35)
        fezeroprime=model.TF_excitatory(fezero+model.external_input_ex_ex, fizero,0.0)
        #fezeroprime=model.TF_excitatory(fezero, fizero,fezero*model.b_e*model.tau_w_e) # with fixe adaptation
        feprimovec[i]=fezeroprime
        fiprimovec[i]=fizero
        # print(i,np.array(xrange1[i]),fezeroprime,fizero)
    import matplotlib.pyplot as plt
    plt.plot(xrange1,feprimovec,'k-',xrange1,xrange1,'k--')
    plt.plot(fiprimovec,feprimovec,'r-')
    plt.show()