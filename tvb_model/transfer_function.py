import numpy as np
import scipy.special as sp_spec

#function from github of Zerlaut :
# https://bitbucket.org/yzerlaut/mean_field_for_multi_input_integration
# ./transfer_functions/theoretical_tools.py
# paper :
#  [ZD_2017] Zerlaut Yann and Destexhe Alain *A mean-field model for conductance-based networks
#        of adaptive exponential integrate-and-fire neurons*, arXiv:1703.00698

#Todo need to add comment on the definition of the function and the parameters of the function like the two last funcstions

class Transfer_Function:
    def __init__(self,model,P_e=None,P_i=None):
        self.model = model
        self.P_e = P_e
        self.P_i = P_i

    def pseq_params(self,P):
        #Todo when P is not define, we can find the good polynome with simulation
        #   for this we can use the function generate_transfer_function in :
        #   https://bitbucket.org/yzerlaut/mean_field_for_multi_input_integration
        #   ./transfer_functions/tf_simulation.py
        if (P is None) or (len(P)==0): # P is empty
            P = [-45e-3]
            for i in range(1,11):
                P.append(0)
                
        # return params['__Q_e'], params['__tau_e'], params['__E_e'], params['__Q_i'], params['__tau_i'], params['__E_i'], params['__g_l'], \
        #        params['__C_m'], params['__E_L'], params['__N_tot'], params['__p_connect'], params['__g'], P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8], P[9], P[10]
        return self.model.Q_e, self.model.tau_e, self.model.E_e, self.model.Q_i, self.model.tau_i, self.model.E_i, self.model.g_l, \
               self.model.C_m, self.model.E_L, self.model.N_tot, self.model.p_connect, self.model.g, P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8], P[9], P[10]


    def get_fluct_regime_vars(self,Fe, Fi, Qe, Te, Ee, Qi, Ti, Ei, Gl, Cm, El, Ntot, pconnec, gei, P0, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10):
        # here TOTAL (sum over synapses) excitatory and inhibitory input

        #firing rate
        fe = Fe*(1.-gei)*pconnec*Ntot # default is 1 !!
        fi = Fi*gei*pconnec*Ntot

        # conductance fluctuation
        # Eqns 5 from [ZD_2017]
        muGe, muGi = Qe*Te*fe, Qi*Ti*fi
        # Eqns 6 from [ZD_2017]
        muG = Gl+muGe+muGi
        muGn, Tm = muG / Gl, Cm / muG

        #membrane potential
        # Eqns 7 from [ZD_2017]
        muV = (muGe*Ee+muGi*Ei+Gl*El)/muG

        # post-synaptic menbrane potential event s around muV
        # Eqns 9 from [ZD_2017]
        Ue, Ui = Qe/muG*(Ee-muV), Qi/muG*(Ei-muV)

        # Standard deviation of the fluctuations
        # Eqns 15 from [ZD_2017]
        sV = np.sqrt(\
                     fe*(Ue*Te)**2/2./(Te+Tm)+\
                     fi*(Qi*Ui)**2/2./(Ti+Tm))

        # Autocorrelation-time of the fluctuations
        # Eqns 17 from [ZD_2017]
        fe, fi = fe+np.finfo(type(fe[0])).eps, fi+np.finfo(type(fi[0])).eps # just to insure a non zero division,
        Tv = ( fe*(Ue*Te)**2 + fi*(Ti*Ui)**2 ) /( fe*(Ue*Te)**2/(Te+Tm) + fi*(Ti*Ui)**2/(Ti+Tm) )
        TvN = Tv*Gl/Cm

        return muV, sV+np.finfo(type(sV[0])).eps, muGn, TvN

    def erfc_func(self,muV, sV, TvN, Vthre, Gl, Cm):
        # Eqns 3 from [ZD_2017]
        return .5/TvN*Gl/Cm*\
          sp_spec.erfc((Vthre-muV)/np.sqrt(2)/sV)

    def threshold_func(self,muV, sV, TvN, muGn, P0, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10):
        """
        setting by default to True the square
        because when use by external modules, coeff[5:]=np.zeros(3)
        in the case of a linear threshold
        """
        # Normalization factors pag. 5 after the equation 4
        muV0, DmuV0 = -60e-3,10e-3
        sV0, DsV0 =4e-3, 6e-3
        TvN0, DTvN0 = 0.5, 1.

        # Eqns 4 from [ZD_2017]
        return P0+P1*(muV-muV0)/DmuV0+\
            P2*(sV-sV0)/DsV0+P3*(TvN-TvN0)/DTvN0+\
            P4*np.log(muGn)+P5*((muV-muV0)/DmuV0)**2+\
            P6*((sV-sV0)/DsV0)**2+P7*((TvN-TvN0)/DTvN0)**2+\
            P8*(muV-muV0)/DmuV0*(sV-sV0)/DsV0+\
            P9*(muV-muV0)/DmuV0*(TvN-TvN0)/DTvN0+\
            P10*(sV-sV0)/DsV0*(TvN-TvN0)/DTvN0

    def TF_my_template(self,fe, fi, Qe, Te, Ee, Qi, Ti, Ei, Gl, Cm, El, Ntot, pconnec, gei, P0, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10):
        # here TOTAL (sum over synapses) excitatory and inhibitory input
        muV, sV, muGn, TvN = self.get_fluct_regime_vars(fe, fi, Qe, Te, Ee, Qi, Ti, Ei, Gl, Cm, El, Ntot, pconnec, gei, P0, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10)
        Vthre = self.threshold_func(muV, sV, TvN, muGn, P0, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10)
        Fout_th = self.erfc_func(muV, sV, TvN, Vthre, Gl, Cm)
        return Fout_th

    def excitatory (self,fe, fi):
        """
        transfer function for excitatory population
        :param fe: firing rate of excitatory population
        :param fi: firing rate of inhibitory population
        :return: result of transfert function
        """
        return self.TF_my_template(fe, fi, *self.pseq_params(self.P_e))

    def inhibitory (self,fe, fi):
        """
        transfer function for inhibitory population
        :param fe: firing rate of excitatory population
        :param fi: firing rate of inhibitory population
        :return: result of transfert function
        """
        return self.TF_my_template(fe, fi, *self.pseq_params(self.P_i))


if __name__=='__main__':

    import numpy as np
    import matplotlib.pylab as plt
    from matplotlib.cm import viridis
    from Zerlaut_model_first_order import Zerlaut_model_first_order
    model = Zerlaut_model_first_order()

    # The model now has the right polynom
    print(model.P_e)

    TF = Transfer_Function(model, model.P_e, model.P_i)

    fig, [axE, axI] = plt.subplots(1, 2, figsize=(8,4))
    fe, fi = np.linspace(0., 20.), np.linspace(0., 20.)
    for i in range(len(fi)):
        axE.plot(fe, TF.excitatory(fe, fi[i]), color=viridis(i/len(fi)))
        axI.plot(fe, TF.inhibitory(fe, fi[i]), color=viridis(i/len(fi)))
    axE.set_title('exc. TF')
    axI.set_title('inh. TF')

    plt.show()
    
