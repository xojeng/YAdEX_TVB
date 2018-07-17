"""
Mean field model based on Master equation about adaptative exponential leacky integrate and fire neurons population
"""

from tvb.simulator.models.base import Model, LOG, numpy, basic, arrays
from transfer_function import Transfer_Function

class Zerlaut_model_second_order(Model):
    r"""
    **References**:
    .. [ZD_2018]  Modeling mesoscopic cortical dynamics using a mean-field model of conductance-based networks of adaptive exponential integrate-and-fire neurons

    Used Eqns 19 from [ZD_2018]_ in ``dfun``.

    #TODO : need to describe a little bit the connectivity between the two populations

    The default parameters are taken from table 1 of [ZD_2018]_, pag.3
    +---------------------------+------------+
    |                 Table 1                |
    +--------------+------------+------------+
    |Parameter     |  Value     | Unit       |
    +==============+============+============+
    |             cellular property          |
    +--------------+------------+------------+
    | g_l          |    1.00    |   nS       |
    +--------------+------------+------------+
    | E_l          |  -65.00    |   mV       |
    +--------------+------------+------------+
    | C_m          |   150.0    |   pF       |
    +--------------+------------+------------+
    | V_thre       |   -50.0    |   mV       |
    +--------------+------------+------------+
    | tau_refrec   |    5.0     |   ms       |
    +--------------+------------+------------+
    | tau_w        |   500.0    |   ms       |
    +--------------+------------+------------+
    |             excitatory cell            |
    +--------------+------------+------------+
    | k_a          |    2.0     |  mV        |
    +--------------+------------+------------+
    | b            |   20.0     |  pA        |
    +--------------+------------+------------+
    | a            |    4.0     |  pA        |
    +--------------+------------+------------+
    |             inhibitory cell            |
    +--------------+------------+------------+
    | k_a          |    0.5     | mV         |
    +--------------+------------+------------+
    | a            |    0.0     | pA         |
    +--------------+------------+------------+
    | b            |    0.0     | ns         |
    +--------------+------------+------------+
    |          synaptic properties           |
    +--------------+------------+------------+
    | E_e          |    0.0     | mV         |
    +--------------+------------+------------+
    | E_i          |   -80.0    | mV         |
    +--------------+------------+------------+
    | Q_e          |    1.0     | nS         |
    +--------------+------------+------------+
    | Q_i          |    5.0     | nS         |
    +--------------+------------+------------+
    | tau_e        |    5.0     | ms         |
    +--------------+------------+------------+
    | tau_i        |    5.0     | ms         |
    +--------------+------------+------------+
    |          numerical network             |
    +--------------+------------+------------+
    | N_tot        |  10000     |            |
    +--------------+------------+------------+
    | epsilon      |    5.0 %   |            |
    +--------------+------------+------------+
    | g            |   20.0 %   |            |
    +--------------+------------+------------+
    | nu_ext       |    4.0     | Hz         |
    +--------------+------------+------------+




    The models (:math:`E`, :math:`I`) phase-plane, including a representation of
    the vector field as well as its nullclines, using default parameters, can be
    seen below:

    .. automethod:: Zerlaut_model_first_order.__init__

    The general formulation for the \textit{\textbf{Zerlaut_model_first_order}} model as a
    dynamical unit at a node $k$ in a BNM with $l$ nodes reads:

    .. math::
         #TODO need to give the equation of the second order

    """
    _ui_name = "Zerlaut_model_2"
    ui_configurable_parameters = ['tau_e', 'tau_i',
                                  'a_e', 'b_e', 'c_e', 'a_i', 'b_i', 'c_i', 'r_e',
                                  'r_i', 'k_e', 'k_i', 'P', 'Q', 'theta_e', 'theta_i',
                                  'alpha_e', 'alpha_i']

    # Define the attributes for this model, these represent possible kwargs.
    g_l = arrays.FloatArray(
        label=":math:`g_{l}`",
        default=numpy.array([10.]), # 10 nS by default, i.e. ~ 100MOhm input resitance at rest
        range=basic.Range(lo=0.1, hi=100.0, step=0.1), # 0.1nS would be a very small cell, 100nS a very big one
        doc="""leak conductance [nS]""",
        order=1)

    E_L = arrays.FloatArray(
        label=":math:`E_{L}`",
        default=numpy.array([-65.0]),
        range=basic.Range(lo=-90.0, hi=-60.0, step=0.1), # resting potential, usually between -85mV and -65mV
        doc="""leak reversal potential [mV]""",
        order=2)

    # N.B. Not independent of g_L, C_m should scale linearly with g_L
    C_m = arrays.FloatArray(
        label=":math:`C_{m}`",
        default=numpy.array([150]),
        range=basic.Range(lo=10.0, hi=500.0, step=10.0), # 20pF very small cell, 400pF very
        doc="""membrane capacitance [pF]""",
        order=3)

    V_th = arrays.FloatArray(
        label=":math:`V_{thre}`",
        default=numpy.array([-50.0]),
        range=basic.Range(lo=-60.0, hi=-40.0, step=0.1),
        doc="""AP threshold [mV]""",
        order=4)

    tau_ref = arrays.FloatArray(
        label=r":math:`\tau_{refrec}`",
        default=numpy.array([5.0]),
        range=basic.Range(lo=0.1, hi=10.0, step=0.1),
        doc="""refractory period [ms]""",
        order=5)

    tau_w = arrays.FloatArray(
        label=r":math:`\tau_w`",
        default=numpy.array([500.0]),
        range=basic.Range(lo=0.0, hi=1000.0, step=10.0), # from 0ms to 1s, sometimes you might want to remove this feature
        doc="""adaptation time constant [ms]""",
        order=6)

    k_a_e = arrays.FloatArray(
        label=":math:`Excitatory k_a`",
        default=numpy.array([2.0]),
        range=basic.Range(lo=0.0, hi=5.0, step=0.1),
        doc="""Excitatory sodium sharpness [mV]""",
        order=7)

    b_e = arrays.FloatArray(
        label=":math:`Excitatory b`",
        default=numpy.array([20.0]),
        range=basic.Range(lo=0.0, hi=100.0, step=0.1), # here also from 0, to remove it if one wants
        doc="""Excitatory adaptation current increment [pA]""",
        order=8)

    a_e = arrays.FloatArray(
        label=":math:`Excitatory a`",
        default=numpy.array([4.0]),
        range=basic.Range(lo=1.0, hi=20.0, step=0.1),
        doc="""Excitatory adaptation conductance [nS]""",
        order=9)

    k_a_i = arrays.FloatArray(
        label=":math:`Inhibitory k_a`",
        default=numpy.array([0.5]),
        range=basic.Range(lo=0.0, hi=5.0, step=0.1),
        doc="""Inhibitory sodium sharpness [mV]""",
        order=10)

    b_i = arrays.FloatArray(
        label=":math:`Inhibitory b`",
        default=numpy.array([0.0]),
        range=basic.Range(lo=0.0, hi=100.0, step=0.1),
        doc="""Inhibitory adaptation current increment [pA]""",
        order=11)

    a_i = arrays.FloatArray(
        label=":math:`Inhibitory a`",
        default=numpy.array([0.0]),
        range=basic.Range(lo=0.0, hi=20.0, step=0.1),
        doc="""Inhibitory adaptation conductance [nS]""",
        order=12)

    E_e = arrays.FloatArray(
        label=r":math:`E_e`",
        default=numpy.array([0.0]),
        range=basic.Range(lo=-20., hi=20., step=0.01),
        doc="""excitatory reversal potential [mV]""",
        order=13)

    E_i = arrays.FloatArray(
        label=":math:`E_i`",
        default=numpy.array([-80.0]),
        range=basic.Range(lo=-100.0, hi=-60.0, step=1.0),
        doc="""inhibitory reversal potential [mV]""",
        order=14)

    Q_e = arrays.FloatArray(
        label=r":math:`Q_e`",
        default=numpy.array([1.0]),
        range=basic.Range(lo=0.0, hi=5.0, step=0.1),
        doc="""excitatory quantal conductance [nS]""",
        order=15)

    Q_i = arrays.FloatArray(
        label=r":math:`Q_i`",
        default=numpy.array([5.0]),
        range=basic.Range(lo=0.0, hi=10.0, step=0.1),
        doc="""inhibitory quantal conductance [nS]""",
        order=16)

    tau_e = arrays.FloatArray(
        label=":math:`\tau_e`",
        default=numpy.array([5.0]),
        range=basic.Range(lo=1.0, hi=10.0, step=1.0),
        doc="""excitatory decay [ms]""",
        order=17)

    tau_i = arrays.FloatArray(
        label=":math:`\tau_i`",
        default=numpy.array([5.0]),
        range=basic.Range(lo=0.5, hi=10.0, step=0.01),
        doc="""inhibitory decay [ms]""",
        order=18)

    N_tot = arrays.IntegerArray(
        label=":math:`N_{tot}`",
        default=numpy.array([10000]),
        range=basic.Range(lo=1000, hi=50000, step=1000),
        doc="""cell number""",
        order=19)

    p_connect = arrays.FloatArray(
        label=":math:`\epsilon`",
        default=numpy.array([0.05]),
        range=basic.Range(lo=0.001, hi=0.2, step=0.001), # valid only for relatively sparse connectivities
        doc="""connectivity probability""",
        order=20)

    g = arrays.FloatArray(
        label=":math:`g`",
        default=numpy.array([0.2]),
        range=basic.Range(lo=0.01, hi=0.4, step=0.01), # inhibitory cell number never overcomes excitatry ones
        doc="""fraction of inhibitory cells""",
        order=21)

    nu_ext = arrays.FloatArray(
        label=":math:`\nu_e^{drive}`",
        default=numpy.array([0.0]),
        range=basic.Range(lo=0.0, hi=100.0, step=0.01),
        doc="""external drive""",
        order=22)

    T = arrays.FloatArray(
        label=":math:`T`",
        default=numpy.array([5.0]),
        range=basic.Range(lo=1.0, hi=20.0, step=0.1),
        doc="""time scale of describing network activity""",
        order=22)

    P_e = arrays.IndexArray(
        label= ":math:`P_e`", #TODO need to check the size of the array when it's used
        # required=False,
        # file_storage=core.FILE_STORAGE_NONE, #TODO need to verified if the value is correct or not
        default= numpy.array([-5.14155819e-02,  6.11046069e-03,  7.44015476e-03,  5.76905434e-05,
                              -1.52902668e-04,  5.63364410e-04,  2.72814525e-04,  5.27870901e-04,
                              -6.77080411e-04,  4.89796440e-04,  1.15073958e-03]),
        doc="""Polynome of excitatory phenomenological threshold (order 10)""",
        order=23)

    P_i = arrays.IndexArray(
        label=":math:`P_i`", #TODO need to check the size of the array when it's used
        # required=False,
        # file_storage=core.FILE_STORAGE_NONE, #TODO need to verified if the value is correct or not
        default= numpy.array([-5.46160664e-02,  4.58415750e-03, -1.77303201e-03,  6.64785219e-04,
                              -3.00637490e-04,  3.93520293e-04, -5.14454957e-04, -6.39186948e-06,
                              -1.39021341e-03, -4.85663596e-04, -3.63617754e-04]),
        doc="""Polynome of inhibitory phenomenological threshold (order 10)""",
        order=24)

    external_input = arrays.FloatArray(
        label=":math:`\nu_e^{drive}`",
        default=numpy.array([0.000]),
        range=basic.Range(lo=0.001, hi=0.1, step=0.001),
        doc="""external drive""",
        order=23)


    # Used for phase-plane axis ranges and to bound random initial() conditions.
    state_variable_range = basic.Dict(
        label="State Variable ranges [lo, hi]",
        default={"E": numpy.array([0.0001, 0.1]), # actually the 200Hz should be replaced by 1/T_refrac, but let's take that
                 "I": numpy.array([0.0001, 0.1]),
                 "C_ee":numpy.array([0.0, 0.0]), # variance is positive or null
                 "C_ei":numpy.array([0.0, 0.0]), # the co-variance is in [-c_ee*c_ii,c_ee*c_ii]
                 "C_ii":numpy.array([0.0, 0.0]) # variance is positive or null
                 },
        doc="""The values for each state-variable should be set to encompass
        the expected dynamic range of that state-variable for the current
        parameters, it is used as a mechanism for bounding random inital
        conditions when the simulation isn't started from an explicit history,
        it is also provides the default range of phase-plane plots.\n
        E: firing rate of excitatory population in KHz\n
        I: firing rate of inhibitory population in KHz\n
        C_ee: the variance of the excitatory population activity \n
        C_ei: the covariance between the excitatory and inhibitory population activities (always symetric) \n
        C_ie: the variance of the inhibitory population activity 
        """,
        order=24)

    variables_of_interest = basic.Enumerate(
        label="Variables watched by Monitors",
        options=["E", "I", "C_ee","C_ei","C_ii"], #TODO need to define the recordable variable
        default=["E"], #TODO need to define the default recordable variables
        select_multiple=True,
        doc="""This represents the default state-variables of this Model to be
               monitored. It can be overridden for each Monitor if desired. The
               corresponding state-variable indices for this model are :math:`E = 0`,
               :math:`I = 1`, :math:`C_ee = 2`, :math:`C_ei = 3`, :math:`C_ii = 4`.""",
        order=25)

    state_variables = 'E I C_ee C_ei C_ii'.split()
    _nvar = 5
    cvar = numpy.array([0,1,2,3,4], dtype=numpy.int32)

    def dfun(self, state_variables, coupling, local_coupling=0.00):
        r"""

        .. math::
           #TODO need to write the equation

        """
        #this variable can be internal varible for optimization
        TF = Transfer_Function(self,self.P_e,self.P_i)
        Ne = self.N_tot * (1-self.g)
        Ni = self.N_tot * self.g

        E = state_variables[0, :]
        I = state_variables[1, :]
        C_ee = state_variables[2, :]
        C_ei = state_variables[3, :]
        C_ii = state_variables[4, :]
        derivative = numpy.empty_like(state_variables)

        # long-range coupling
        c_0 = coupling[0, :]

        # short-range (local) coupling
        lc_E = local_coupling * E +self.external_input
        lc_I = local_coupling * I

        #equation is inspired from github of Zerlaut :
        # function from github of Zerlaut :
        # https://github.com/yzerlaut/notebook_papers/blob/master/modeling_mesoscopic_dynamics/mean_field/master_equation.py

        # Excitatory firing rate derivation
        # #TODO : need to check the accuracy of all equations
        derivative[0] = 1./self.T*(\
                .5*C_ee*self._diff2_fe_fe(TF.excitatory, E+c_0+lc_E, I+lc_I)+\
                .5*C_ei*self._diff2_fe_fi(TF.excitatory, E+c_0+lc_E, I+lc_I)+\
                .5*C_ei*self._diff2_fi_fe(TF.excitatory, E+c_0+lc_E, I+lc_I)+\
                .5*C_ii*self._diff2_fi_fi(TF.excitatory, E+c_0+lc_E, I+lc_I)+\
                TF.excitatory(E+c_0+lc_E, I+lc_I)-E)
        #Inhibitory firing rate derivation
        derivative[1] = 1./self.T*(\
                .5*C_ee*self._diff2_fe_fe(TF.inhibitory, E, I+lc_I)+\
                .5*C_ei*self._diff2_fe_fi(TF.inhibitory, E, I+lc_I)+\
                .5*C_ei*self._diff2_fi_fe(TF.inhibitory, E, I+lc_I)+\
                .5*C_ii*self._diff2_fi_fi(TF.inhibitory, E, I+lc_I)+\
                TF.inhibitory(E, I+lc_I)-I)
        #Covariance excitatory-excitatory derivation
        derivative[2] = 1./self.T*(\
                1./Ne*TF.excitatory(E+c_0+lc_E, I+lc_I)*\
                (1./self.T-TF.excitatory(E+c_0+lc_E, I+lc_I))+\
                (TF.excitatory(E+c_0+lc_E, I+lc_I)-E)**2+\
                2.*C_ee*self._diff_fe(TF.excitatory, E+c_0+lc_E, I+lc_I)+\
                2.*C_ei*self._diff_fi(TF.inhibitory, E, I+lc_I)+\
                -2.*C_ee)


        #Covariance excitatory-inhibitory or inhibitory-excitatory derivation
        derivative[3] = 1./self.T*\
                ((TF.excitatory(E+c_0+lc_E, I+lc_I)-E)*\
                (TF.inhibitory(E+c_0+lc_E, I+lc_I)-I)+\
                C_ee*self._diff_fe(TF.excitatory, E+c_0+lc_E, I+lc_I)+\
                C_ei*self._diff_fe(TF.inhibitory, E, I+lc_I)+\
                C_ei*self._diff_fi(TF.excitatory, E+c_0+lc_E, I+lc_I)+\
                C_ii*self._diff_fi(TF.inhibitory, E, I+lc_I)+\
                -2.*C_ei)

        #Covariance inhibitory-inhibitory derivation
        derivative[4] = 1./self.T*(\
                1./Ni*TF.inhibitory(E, I+lc_I)*\
                (1./self.T-TF.inhibitory(E, I+lc_I))+\
                (TF.inhibitory(E, I+lc_I)-I)**2+\
                2.*C_ii*self._diff_fi(TF.inhibitory, E, I+lc_I)+\
                2.*C_ei*self._diff_fe(TF.excitatory, E+c_0+lc_E, I+lc_I)+\
                -2.*C_ii)

        return derivative

    #equation is inspired from github of Zerlaut :
    # function from github of Zerlaut :
    # https://github.com/yzerlaut/notebook_papers/blob/master/modeling_mesoscopic_dynamics/mean_field/master_equation.py
    ##### Derivatives taken numerically : use a central difference formula with spacing `dx`
    ## to be implemented analitically ! not hard...

    def _diff_fe(self,TF, fe, fi,df=1e-7):
        return (TF(fe+df, fi)-TF(fe-df,fi))/(2*df*1e3)

    def _diff_fi(self,TF, fe, fi,df=1e-7):
        return (TF(fe, fi+df)-TF(fe, fi-df))/(2*df*1e3)

    def _diff2_fe_fe(self,TF, fe, fi,df=1e-7):
        return (TF(fe+df, fi)-2*TF(fe,fi)+TF(fe-df, fi))/((df*1e3)**2)

    def _diff2_fi_fe(self,TF, fe, fi,df=1e-7):
        return (self._diff_fi(TF, fe+df, fi)-self._diff_fi(TF,fe-df,fi))/(2*df*1e3)

    def _diff2_fe_fi(self,TF, fe, fi,df=1e-7):
        return (self._diff_fe(TF, fe, fi+df)-self._diff_fe(TF,fe, fi-df))/(2*df*1e3)

    def _diff2_fi_fi(self,TF, fe, fi,df=1e-7):
        return (TF(fe, fi+df)-2*TF(fe, fi)+TF(fe, fi-df))/((df*1e3)**2)
