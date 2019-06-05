from tvb.simulator.lab import *
import numpy as np
import matplotlib.pyplot as plt
from tvb.datatypes.time_series import TimeSeriesRegion
import tvb.analyzers.correlation_coefficient as corr_coeff
import tvb.analyzers.fcd_matrix as fcd
import matplotlib
matplotlib.rcParams['figure.figsize'] = (20.0, 10.0)
import tools_simulation as tools


def plot_bif_diagram(rangeG, parameter_simulation, plot_name=None, save=False):

    result = {'low':[], 'high':[]}
    for G in rangeG:
        for initconds in ['low', 'high']:
            if initconds == 'low':
                path = parameter_simulation['path_result']+"low_G_%s/"%str(G)
            elif initconds == 'high':
                path = parameter_simulation['path_result']+"high_G_%s/"%str(G)

            result[initconds].append(tools.get_result(path, 4999., 5000.)[0][1][-1])

    Fe_low = np.array(result['low'])[:,0].max(axis=1)
    Fe_high = np.array(result['high'])[:,0].max(axis=1)
    Fi_low = np.array(result['low'])[:,1].max(axis=1)
    Fi_high = np.array(result['high'])[:,1].max(axis=1)

    plt.figure(figsize=(8, 8))
    # plot low activity
    plt.plot(1e3*Fe_low, 'o', markeredgecolor='b', markerfacecolor='none', markersize=10, label='low act. exc')
    plt.plot(1e3*Fi_low, 'o', markeredgecolor='c', markerfacecolor='none', markersize=10, label='low act. inh')

    # plot high activity
    plt.plot(1e3*Fe_high, 'ro', markeredgecolor='none', markersize=6, label='high act. exc')
    plt.plot(1e3*Fi_high, 'mo', markeredgecolor='none', markersize=6, label='high act. inh')

    plt.xticks(np.arange(len(rangeG)), rangeG)
    plt.xlabel('$G_{coupl}$', fontsize=20)
    plt.ylabel('Firing rate (in Hz)', fontsize=20)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    if plot_name == None:
        title = 'bifurcation_diagram_zerlaut'
    else:
        title = 'bifurcation_diagram_zerlaut_'+plot_name

    plt.title(title, fontsize=20)
    plt.tight_layout()
    if save == True:
        plt.savefig('/home/mansi/summer_project_2019/zerlaut/figures/%s.png'%(title))
    else:
        pass

    plt.show()