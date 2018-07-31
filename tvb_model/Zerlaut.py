"""
Mean field model based on Master equation about adaptative exponential leacky integrate and fire neurons population
"""

from tvb.simulator.models.base import Model, LOG, numpy, basic, arrays, core
import scipy.special as sp_spec


class Zerlaut_first_order(Model):
    r"""
    **References**:
    .. [ZD_2018]  Modeling mesoscopic cortical dynamics using a mean-field model of conductance-based networks of adaptive exponential integrate-and-fire neurons

    Used Eqns 20 from [ZD_2018]_ in ``dfun``.

    The default parameters are taken from table 1 of [ZD_20108]_, pag.47
    +---------------------------+------------+
    |                 Table 1                |
    +--------------+------------+------------+
    |Parameter     |  Value     | Unit       |
    +==============+============+============+
    |             cellular property          |
    +--------------+------------+------------+
    | g_L          |    1.00    |   nS       |
    +--------------+------------+------------+
    | E_L          |  -65.00    |   mV       |
    +--------------+------------+------------+
    | C_m          |   150.0    |   pF       |
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
    | p_connect    |    5.0 %   |            |
    +--------------+------------+------------+
    | g            |   20.0 %   |            |
    +--------------+------------+------------+
    |external_input|    0.004   | Hz         |
    +--------------+------------+------------+

    The default coefficients of the transfert function are taken from table 2 of [ZD_20108]_, pag.49
    +---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+
    |             excitatory cell   |
    +---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+
    |-5.14155819e-02| 6.11046069e-03| 7.44015476e-03| 5.76905434e-05|-1.52902668e-04| 5.63364410e-04| 2.72814525e-04| 5.27870901e-04|-6.77080411e-04| 4.89796440e-04| 1.15073958e-03|
    +---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+
    |             inhibitory cell   |
    +---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+
    |-5.46160664e-02| 4.58415750e-03|-1.77303201e-03| 6.64785219e-04|-3.00637490e-04| 3.93520293e-04|-5.14454957e-04|-6.39186948e-06|-1.39021341e-03|-4.85663596e-04|-3.63617754e-04|
    +---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+---------------+

    The models (:math:`E`, :math:`I`) phase-plane, including a representation of
    the vector field as well as its nullclines, using default parameters, can be
    seen below:

    .. automethod:: Zerlaut_first_order.__init__

    The general formulation for the \textit{\textbf{Zerlaut_first_order}} model as a
    dynamical unit at a node $k$ in a BNM with $l$ nodes reads:

    .. math::
            T\dot{E}_k &= F_e-E_k  \\
            T\dot{I}_k &= F_i-I_k  \\
            F_\lambda = Erfc(V^{eff}_{thre}-\mu_V/\sqrt(2)\sigma_V)

    """
    _ui_name = "Zerlaut_first_order"
    ui_configurable_parameters = ['g_L', 'E_L', 'C_m',
                                  'E_e', 'E_i', 'Q_e', 'Q_i', 'tau_e', 'tau_i',
                                  'N_tot', 'p_connect', 'g', 'T',
                                  'external_input']

    # Define traited attributes for this model, these represent possible kwargs.
    g_L = arrays.FloatArray(
        label=":math:`g_{L}`",
        default=numpy.array([10.]),  # 10 nS by default, i.e. ~ 100MOhm input resitance at rest
        range=basic.Range(lo=0.1, hi=100.0, step=0.1),  # 0.1nS would be a very small cell, 100nS a very big one
        doc="""leak conductance [nS]""",
        order=1)

    E_L = arrays.FloatArray(
        label=":math:`E_{L}`",
        default=numpy.array([-65.0]),
        range=basic.Range(lo=-90.0, hi=-60.0, step=0.1),  # resting potential, usually between -85mV and -65mV
        doc="""leak reversal potential [mV]""",
        order=2)

    # N.B. Not independent of g_L, C_m should scale linearly with g_L
    C_m = arrays.FloatArray(
        label=":math:`C_{m}`",
        default=numpy.array([150]),
        range=basic.Range(lo=10.0, hi=500.0, step=10.0),  # 20pF very small cell, 400pF very
        doc="""membrane capacitance [pF]""",
        order=3)

    E_e = arrays.FloatArray(
        label=r":math:`E_e`",
        default=numpy.array([0.0]),
        range=basic.Range(lo=-20., hi=20., step=0.01),
        doc="""excitatory reversal potential [mV]""",
        order=4)

    E_i = arrays.FloatArray(
        label=":math:`E_i`",
        default=numpy.array([-80.0]),
        range=basic.Range(lo=-100.0, hi=-60.0, step=1.0),
        doc="""inhibitory reversal potential [mV]""",
        order=5)

    Q_e = arrays.FloatArray(
        label=r":math:`Q_e`",
        default=numpy.array([1.0]),
        range=basic.Range(lo=0.0, hi=5.0, step=0.1),
        doc="""excitatory quantal conductance [nS]""",
        order=6)

    Q_i = arrays.FloatArray(
        label=r":math:`Q_i`",
        default=numpy.array([5.0]),
        range=basic.Range(lo=0.0, hi=10.0, step=0.1),
        doc="""inhibitory quantal conductance [nS]""",
        order=7)

    tau_e = arrays.FloatArray(
        label=":math:`\tau_e`",
        default=numpy.array([5.0]),
        range=basic.Range(lo=1.0, hi=10.0, step=1.0),
        doc="""excitatory decay [ms]""",
        order=8)

    tau_i = arrays.FloatArray(
        label=":math:`\tau_i`",
        default=numpy.array([5.0]),
        range=basic.Range(lo=0.5, hi=10.0, step=0.01),
        doc="""inhibitory decay [ms]""",
        order=9)

    N_tot = arrays.IntegerArray(
        label=":math:`N_{tot}`",
        default=numpy.array([10000]),
        range=basic.Range(lo=1000, hi=50000, step=1000),
        doc="""cell number""",
        order=10)

    p_connect = arrays.FloatArray(
        label=":math:`\epsilon`",
        default=numpy.array([0.05]),
        range=basic.Range(lo=0.001, hi=0.2, step=0.001),  # valid only for relatively sparse connectivities
        doc="""connectivity probability""",
        order=11)

    g = arrays.FloatArray(
        label=":math:`g`",
        default=numpy.array([0.2]),
        range=basic.Range(lo=0.01, hi=0.4, step=0.01),  # inhibitory cell number never overcomes excitatry ones
        doc="""fraction of inhibitory cells""",
        order=12)

    T = arrays.FloatArray(
        label=":math:`T`",
        default=numpy.array([5.0]),
        range=basic.Range(lo=1., hi=20.0, step=0.1),
        doc="""time scale of describing network activity""",
        order=13)

    P_e = arrays.IndexArray(
        label=":math:`P_e`",  # TODO need to check the size of the array when it's used
        default=numpy.array([-5.14155819e-02,  6.11046069e-03,  7.44015476e-03,  5.76905434e-05,
                              -1.52902668e-04,  5.63364410e-04,  2.72814525e-04,  5.27870901e-04,
                              -6.77080411e-04,  4.89796440e-04,  1.15073958e-03]),
        doc="""Polynome of excitatory phenomenological threshold (order 10)""",
        order=14)

    P_i = arrays.IndexArray(
        label=":math:`P_i`",  # TODO need to check the size of the array when it's used
        default=numpy.array([-5.46160664e-02,  4.58415750e-03, -1.77303201e-03,  6.64785219e-04,
                              -3.00637490e-04,  3.93520293e-04, -5.14454957e-04, -6.39186948e-06,
                              -1.39021341e-03, -4.85663596e-04, -3.63617754e-04]),
        doc="""Polynome of inhibitory phenomenological threshold (order 10)""",
        order=15)

    external_input = arrays.FloatArray(
        label=":math:`\nu_e^{drive}`",
        default=numpy.array([0.004]),
        range=basic.Range(lo=0.001, hi=0.1, step=0.001),
        doc="""external drive""",
        order=16)

    # Used for phase-plane axis ranges and to bound random initial() conditions.
    state_variable_range = basic.Dict(
        label="State Variable ranges [lo, hi]",
        default={"E": numpy.array([0.0, 0.2]),  # actually the 200Hz should be replaced by 1/T_refrac
                 "I": numpy.array([0.0, 0.2])},
        doc="""The values for each state-variable should be set to encompass
        the expected dynamic range of that state-variable for the current
        parameters, it is used as a mechanism for bounding random initial
        conditions when the simulation isn't started from an explicit history,
        it is also provides the default range of phase-plane plots.\n
        E: firing rate of excitatory population in KHz\n
        I: firing rate of inhibitory population in kHz
        """,
        order=17)

    variables_of_interest = basic.Enumerate(
        label="Variables watched by Monitors",
        options=["E", "I"],
        default=["E"],
        select_multiple=True,
        doc="""This represents the default state-variables of this Model to be
               monitored. It can be overridden for each Monitor if desired. The
               corresponding state-variable indices for this model are :math:`E = 0`
               and :math:`I = 1`.""",
        order=18)

    state_variables = 'E I'.split()
    _nvar = 2
    cvar = numpy.array([0, 1], dtype=numpy.int32)

    def dfun(self, state_variables, coupling, local_coupling=0.00):
        r"""
        .. math::
            T \dot{\nu_\mu} &= -F_\mu(\nu_e,\nu_i) + \nu_\mu ,\all\mu\in\{e,i\}

        """
        E = state_variables[0, :]
        I = state_variables[1, :]
        derivative = numpy.empty_like(state_variables)

        # long-range coupling
        c_0 = coupling[0, :]

        # short-range (local) coupling
        lc_E = local_coupling * E
        lc_I = local_coupling * I

        # Excitatory firing rate derivation
        derivative[0] = (self.TF_excitatory(E+c_0+lc_E+self.external_input, I+lc_I+self.external_input)-E)/self.T
        # Inhibitory firing rate derivation
        derivative[1] = (self.TF_inhibitory(E+self.external_input, I+lc_I+self.external_input)-I)/self.T

        return derivative

    def TF_excitatory(self, fe, fi):
        """
        transfer function for excitatory population
        :param fe: firing rate of excitatory population
        :param fi: firing rate of inhibitory population
        :return: result of transfer function
        """
        return self.TF(fe, fi, self.P_e)

    def TF_inhibitory(self, fe, fi):
        """
        transfer function for inhibitory population
        :param fe: firing rate of excitatory population
        :param fi: firing rate of inhibitory population
        :return: result of transfer function
        """
        return self.TF(fe, fi, self.P_i)

    def TF(self, fe, fi, P):
        """
        transfer function for inhibitory population
        Inspired from the next repository :
        https://github.com/yzerlaut/notebook_papers/tree/master/modeling_mesoscopic_dynamics
        :param fe: firing rate of excitatory population
        :param fi: firing rate of inhibitory population
        :param P: Polynome of neurons phenomenological threshold (order 10)
        :return: result of transfer function
        """
        mu_V, sigma_V, mu_Gn, T_V = self.get_fluct_regime_vars(fe, fi, self.Q_e, self.tau_e, self.E_e,
                                                               self.Q_i, self.tau_i, self.E_i,
                                                               self.g_L, self.C_m, self.E_L, self.N_tot,
                                                               self.p_connect, self.g)
        V_thre = self.threshold_func(mu_V, sigma_V, T_V*self.g_L/self.C_m, mu_Gn,
                                     P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8], P[9], P[10])
        V_thre *= 1e3  # the threshold need to be in mv and not in Volt
        f_out = self.estimate_firing_rate(mu_V, sigma_V, T_V, V_thre)
        return f_out

    @staticmethod
    def get_fluct_regime_vars(Fe, Fi, Q_e, tau_e, E_e, Q_i, tau_i, E_i, g_L, C_m, E_L, N_tot, p_connect, g):
        # firing rate
        # 1e-6 represent spontaneous release of synaptic neurotransmitter or some intrinsic currents of neurons
        fe = (Fe+1e-6)*(1.-g)*p_connect*N_tot
        fi = (Fi+1e-6)*g*p_connect*N_tot
        # conductance fluctuation and effective membrane time constant
        mu_Ge, mu_Gi = Q_e*tau_e*fe, Q_i*tau_i*fi  # Eqns 5 from [ZD_2018]
        mu_G = g_L+mu_Ge+mu_Gi  # Eqns 6 from [ZD_2018]
        mu_Gn, T_m = mu_G/g_L, C_m/mu_G
        # membrane potential
        mu_V = (mu_Ge*E_e+mu_Gi*E_i+g_L*E_L)/mu_G  # Eqns 7 from [ZD_2018]
        # post-synaptic membrane potential event s around muV
        U_e, U_i = Q_e/mu_G*(E_e-mu_V), Q_i/mu_G*(E_i-mu_V)  # Eqns 9 from [ZD_2018]
        # Standard deviation of the fluctuations
        # Eqns 15 from [ZD_2018]
        sigma_V = numpy.sqrt(fe*(U_e*tau_e)**2/(2.*(tau_e+T_m))+fi*(U_i*tau_i)**2/(2.*(tau_i+T_m)))
        # Autocorrelation-time of the fluctuations Eqns 17 from [ZD_2018]
        T_V_numerator = (fe*(U_e*tau_e)**2 + fi*(U_i*tau_i)**2)
        T_V_denominator = (fe*(U_e*tau_e)**2/(tau_e+T_m) + fi*(U_i*tau_i)**2/(tau_i+T_m))
        T_V = numpy.divide(T_V_numerator, T_V_denominator, out=numpy.ones_like(T_V_numerator),
                           where=T_V_denominator != 0.0)
        return mu_V, sigma_V, mu_Gn, T_V

    @staticmethod
    def threshold_func(muV, sigmaV, TvN, muGn, P0, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10):
        # Normalization factors page 48 after the equation 4 from [ZD_2018]
        muV0, DmuV0 = -60.0, 10.0
        sV0, DsV0 = 4.0, 6.0
        TvN0, DTvN0 = 0.5, 1.
        V = (muV-muV0)/DmuV0
        S = (sigmaV-sV0)/DsV0
        T = (TvN-TvN0)/DTvN0
        # Eqns 4 from [ZD_2018]
        return P0 + P1*V + P2*S + P3*T + P4*numpy.log(muGn) + P5*V**2 + P6*S**2 + P7*T**2 + P8*V*S + P9*V*T + P10*S*T

    @staticmethod
    def estimate_firing_rate(muV, sigmaV, Tv, Vthre):
        # Eqns 3 from [ZD_2018]
        return sp_spec.erfc((Vthre-muV) / (numpy.sqrt(2)*sigmaV)) / (2*Tv)


class Zerlaut_second_order(Zerlaut_first_order):
    r"""
    **References**:
    .. [ZD_2018]  Modeling mesoscopic cortical dynamics using a mean-field model of conductance-based networks of adaptive exponential integrate-and-fire neurons

    Used Eqns 18 from [ZD_2018]_ in ``dfun``.

    (See Zerlaut_first_order for the default value)

    The models (:math:`E`, :math:`I`) phase-plane, including a representation of
    the vector field as well as its nullclines, using default parameters, can be
    seen below:

    .. automethod:: Zerlaut_second_order.__init__

    The general formulation for the \textit{\textbf{Zerlaut_second_order}} model as a
    dynamical unit at a node $k$ in a BNM with $l$ nodes reads:

    .. math::
        \forall \mu,\lambda,\eta \in \{e,i\}^3\, ,
        \left\{
        \begin{split}
        T \, \frac{\partial \nu_\mu}{\partial t} = & (\mathcal{F}_\mu - \nu_\mu )
        + \frac{1}{2} \, c_{\lambda \eta} \,
        \frac{\partial^2 \mathcal{F}_\mu}{\partial \nu_\lambda \partial \nu_\eta} \\
        T \, \frac{\partial c_{\lambda \eta} }{\partial t}  =  & A_{\lambda \eta} +
        (\mathcal{F}_\lambda - \nu_\lambda ) \, (\mathcal{F}_\eta - \nu_\eta ) + \\
        & c_{\lambda \mu} \frac{\partial \mathcal{F}_\mu}{\partial \nu_\lambda} +
        c_{\mu \eta} \frac{\partial \mathcal{F}_\mu}{\partial \nu_\eta}
        - 2  c_{\lambda \eta}
        \end{split}
        \right.

        with:
        A_{\lambda \eta} =
        \left\{
        \begin{split}
        \frac{\mathcal{F}_\lambda \, (1/T - \mathcal{F}_\lambda)}{N_\lambda}
        \qquad & \textrm{if  } \lambda=\eta \\
        0 \qquad & \textrm{otherwise}
        \end{split}
        \right.
    """

    _ui_name = "Zerlaut_second_order"

    #  Used for phase-plane axis ranges and to bound random initial() conditions.
    state_variable_range = basic.Dict(
        label="State Variable ranges [lo, hi]",
        default={"E": numpy.array([0.0, 0.2]), # actually the 200Hz should be replaced by 1/T_refrac
                 "I": numpy.array([0.0, 0.2]),
                 "C_ee": numpy.array([0.0, 0.0]),  # variance is positive or null
                 "C_ei": numpy.array([0.0, 0.0]),  # the co-variance is in [-c_ee*c_ii,c_ee*c_ii]
                 "C_ii": numpy.array([0.0, 0.0])  # variance is positive or null
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
        order=17)

    variables_of_interest = basic.Enumerate(
        label="Variables watched by Monitors",
        options=["E", "I", "C_ee","C_ei","C_ii"],
        default=["E"],
        select_multiple=True,
        doc="""This represents the default state-variables of this Model to be
               monitored. It can be overridden for each Monitor if desired. The
               corresponding state-variable indices for this model are :math:`E = 0`,
               :math:`I = 1`, :math:`C_ee = 2`, :math:`C_ei = 3`, :math:`C_ii = 4`.""",
        order=18)

    state_variables = 'E I C_ee C_ei C_ii'.split()
    _nvar = 5
    cvar = numpy.array([0, 1, 2, 3, 4], dtype=numpy.int32)

    def dfun(self, state_variables, coupling, local_coupling=0.00):
        r"""
        .. math::
            \forall \mu,\lambda,\eta \in \{e,i\}^3\, ,
            \left\{
            \begin{split}
            T \, \frac{\partial \nu_\mu}{\partial t} = & (\mathcal{F}_\mu - \nu_\mu )
            + \frac{1}{2} \, c_{\lambda \eta} \,
            \frac{\partial^2 \mathcal{F}_\mu}{\partial \nu_\lambda \partial \nu_\eta} \\
            T \, \frac{\partial c_{\lambda \eta} }{\partial t}  =  & A_{\lambda \eta} +
            (\mathcal{F}_\lambda - \nu_\lambda ) \, (\mathcal{F}_\eta - \nu_\eta ) + \\
            & c_{\lambda \mu} \frac{\partial \mathcal{F}_\mu}{\partial \nu_\lambda} +
            c_{\mu \eta} \frac{\partial \mathcal{F}_\mu}{\partial \nu_\eta}
            - 2  c_{\lambda \eta}
            \end{split}
            \right.

            with:
            A_{\lambda \eta} =
            \left\{
            \begin{split}
            \frac{\mathcal{F}_\lambda \, (1/T - \mathcal{F}_\lambda)}{N_\lambda}
            \qquad & \textrm{if  } \lambda=\eta \\
            0 \qquad & \textrm{otherwise}
            \end{split}
            \right.

        """
        # this variable can be internal variable for optimization
        N_e = self.N_tot * (1-self.g)
        N_i = self.N_tot * self.g

        E = state_variables[0, :]
        I = state_variables[1, :]
        C_ee = state_variables[2, :]
        C_ei = state_variables[3, :]
        C_ii = state_variables[4, :]
        derivative = numpy.empty_like(state_variables)

        # long-range coupling
        c_0 = coupling[0, :]

        # short-range (local) coupling
        lc_E = local_coupling * E
        lc_I = local_coupling * I

        E_input_excitatory = E+c_0+lc_E+self.external_input
        E_input_inhibitory = E+lc_E+self.external_input
        I_input_excitatory = I+lc_I+self.external_input
        I_input_inhibitory = I+lc_I+self.external_input

        _TF_e = self.TF_excitatory(E_input_excitatory, I_input_excitatory)
        _TF_i = self.TF_inhibitory(E_input_inhibitory, I_input_inhibitory)

        # Derivatives taken numerically : use a central difference formula with spacing `dx`
        def _diff_fe(TF, fe, fi, df=1e-7):
            return (TF(fe+df, fi)-TF(fe-df, fi))/(2*df*1e3)

        def _diff_fi(TF, fe, fi, df=1e-7):
            return (TF(fe, fi+df)-TF(fe, fi-df))/(2*df*1e3)

        def _diff2_fe_fe_e(fe, fi, df=1e-7):
            TF = self.TF_excitatory
            return (TF(fe+df, fi)-2*_TF_e+TF(fe-df, fi))/((df*1e3)**2)

        def _diff2_fe_fe_i(fe, fi, df=1e-7):
            TF = self.TF_inhibitory
            return (TF(fe+df, fi)-2*_TF_i+TF(fe-df, fi))/((df*1e3)**2)

        def _diff2_fi_fe(TF, fe, fi, df=1e-7):
            return (_diff_fi(TF, fe+df, fi)-_diff_fi(TF, fe-df, fi))/(2*df*1e3)

        def _diff2_fe_fi(TF, fe, fi, df=1e-7):
            return (_diff_fe(TF, fe, fi+df)-_diff_fe(TF, fe, fi-df))/(2*df*1e3)

        def _diff2_fi_fi_e(fe, fi, df=1e-7):
            TF = self.TF_excitatory
            return (TF(fe, fi+df)-2*_TF_e+TF(fe, fi-df))/((df*1e3)**2)

        def _diff2_fi_fi_i(fe, fi, df=1e-7):
            TF = self.TF_inhibitory
            return (TF(fe, fi+df)-2*_TF_i+TF(fe, fi-df))/((df*1e3)**2)

        _diff_fe_TF_e = _diff_fe(self.TF_excitatory, E_input_excitatory, I_input_excitatory)
        _diff_fe_TF_i = _diff_fe(self.TF_inhibitory, E_input_inhibitory, I_input_inhibitory)
        _diff_fi_TF_e = _diff_fi(self.TF_excitatory, E_input_excitatory, I_input_excitatory)
        _diff_fi_TF_i = _diff_fi(self.TF_inhibitory, E_input_inhibitory, I_input_inhibitory)

        # equation is inspired from github of Zerlaut :
        # https://github.com/yzerlaut/notebook_papers/blob/master/modeling_mesoscopic_dynamics/mean_field/master_equation.py
        # Excitatory firing rate derivation
        derivative[0] = (_TF_e - E
                         + .5*C_ee*_diff2_fe_fe_e(E_input_excitatory, I_input_excitatory)
                         + .5*C_ei*_diff2_fe_fi(self.TF_excitatory, E_input_excitatory, I_input_excitatory)
                         + .5*C_ei*_diff2_fi_fe(self.TF_excitatory, E_input_excitatory, I_input_excitatory)
                         + .5*C_ii*_diff2_fi_fi_e(E_input_excitatory, I_input_excitatory)
                         )/self.T
        # Inhibitory firing rate derivation
        derivative[1] = (_TF_i - I
                         + .5*C_ee*_diff2_fe_fe_i(E_input_inhibitory, I_input_inhibitory)
                         + .5*C_ei*_diff2_fe_fi(self.TF_inhibitory, E_input_inhibitory, I_input_inhibitory)
                         + .5*C_ei*_diff2_fi_fe(self.TF_inhibitory, E_input_inhibitory, I_input_inhibitory)
                         + .5*C_ii*_diff2_fi_fi_i(E_input_inhibitory, I_input_inhibitory)
                         )/self.T
        # Covariance excitatory-excitatory derivation
        derivative[2] = (_TF_e*(1./self.T-_TF_e)/N_e
                         + (_TF_e-E)**2
                         + 2.*C_ee*_diff_fe_TF_e
                         + 2.*C_ei*_diff_fi_TF_i
                         - 2.*C_ee
                         )/self.T
        # Covariance excitatory-inhibitory or inhibitory-excitatory derivation
        derivative[3] = ((_TF_e-E)*(_TF_i-I)
                         + C_ee*_diff_fe_TF_e
                         + C_ei*_diff_fe_TF_i
                         + C_ei*_diff_fi_TF_e
                         + C_ii*_diff_fi_TF_i
                         - 2.*C_ei
                         )/self.T
        # Covariance inhibitory-inhibitory derivation
        derivative[4] = (_TF_i*(1./self.T-_TF_i)/N_i
                         + (_TF_i-I)**2
                         + 2.*C_ii*_diff_fi_TF_i
                         + 2.*C_ei*_diff_fe_TF_e
                         - 2.*C_ii
                         )/self.T

        return derivative

