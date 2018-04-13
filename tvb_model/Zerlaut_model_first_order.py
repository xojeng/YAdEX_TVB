"""
Mean field model based on Master equation about adaptative exponential leacky integrate and fire neurons population
"""

from tvb.simulator.models.base import Model, LOG, numpy, basic, arrays, core
from transfer_function import Transfer_Function

class Zerlaut_model_first_order(Model):
    r"""
    **References**:
    .. [ZD_2017] Zerlaut Yann and Destexhe Alain *A mean-field model for conductance-based networks
        of adaptive exponential integrate-and-fire neurons*, arXiv:1703.00698

    Used Eqns 20 from [ZD_2017]_ in ``dfun``.

    #TODO : need to discribe a little bit the connectivity between two populations

    The default parameters are taken from table 1 of [ZD_20107]_, pag.3
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
            T\dot{E}_k &= F_e-E_k  \\
            T\dot{I}_k &= F_i-I_k  \\
            F_\lambda = #TODO need to give the expression of F

    """
    _ui_name = "Zerlaut_model_1"
    ui_configurable_parameters = ['c_ee', 'c_ei', 'c_ie', 'c_ii', 'tau_e', 'tau_i',
                                  'a_e', 'b_e', 'c_e', 'a_i', 'b_i', 'c_i', 'r_e',
                                  'r_i', 'k_e', 'k_i', 'P', 'Q', 'theta_e', 'theta_i',
                                  'alpha_e', 'alpha_i']

    # Define traited attributes for this model, these represent possible kwargs.
    g_l = arrays.FloatArray(
        label=":math:`g_{l}`",
        default=numpy.array([1.0]),
        range=basic.Range(lo=0.0, hi=10.0, step=0.1), #TODO need to define the range
        doc="""leak conductance [nS]""",
        order=1)

    E_L = arrays.FloatArray(
        label=":math:`E_{L}`",
        default=numpy.array([-65.0]),
        range=basic.Range(lo=-70.0, hi=-60.0, step=0.1), #TODO need to define the range
        doc="""leak reversal potential [mV]""",
        order=2)

    C_m = arrays.FloatArray(
        label=":math:`C_{m}`",
        default=numpy.array([150]),
        range=basic.Range(lo=100.0, hi=200.0, step=10.0), #TODO need to define the range
        doc="""membrane capacitance [pF]""",
        order=3)

    V_th = arrays.FloatArray(
        label=":math:`V_{thre}`",
        default=numpy.array([-50.0]),
        range=basic.Range(lo=-60.0, hi=-40.0, step=0.1), #TODO need to define the range
        doc="""AP threshold [mV]""",
        order=4)

    tau_ref = arrays.FloatArray(
        label=r":math:`\tau_{refrec}`",
        default=numpy.array([5.0]),
        range=basic.Range(lo=0.0, hi=10.0, step=0.01), #TODO need to define the range
        doc="""refractory period [ms]""",
        order=5)

    tau_w = arrays.FloatArray(
        label=r":math:`\tau_w`",
        default=numpy.array([500.0]),
        range=basic.Range(lo=400.0, hi=600.0, step=10.0), #TODO need to define the range
        doc="""adaptation time constant [ms]""",
        order=6)

    k_a_e = arrays.FloatArray(
        label=":math:`Excitatory k_a`",
        default=numpy.array([2.0]),
        range=basic.Range(lo=0.0, hi=3.0, step=0.1), #TODO need to define the range
        doc="""Excitatory sodium sharpness [mV]""",
        order=7)

    b_e = arrays.FloatArray(
        label=":math:`Excitatory b`",
        default=numpy.array([20.0]),
        range=basic.Range(lo=10.0, hi=30.0, step=0.1), #TODO need to define the range
        doc="""Excitatory adaptation current increment [pA]""",
        order=8)

    a_e = arrays.FloatArray(
        label=":math:`Excitatory a`",
        default=numpy.array([4.0]),
        range=basic.Range(lo=1.0, hi=20.0, step=0.1), #TODO need to define the range
        doc="""Excitatory adaptation conductance [nS]""",
        order=9)

    k_a_i = arrays.FloatArray(
        label=":math:`Inhibitory k_a`",
        default=numpy.array([0.5]),
        range=basic.Range(lo=0.0, hi=3.0, step=0.1), #TODO need to define the range
        doc="""Inhibitory sodium sharpness [mV]""",
        order=10)

    b_i = arrays.FloatArray(
        label=":math:`Inhibitory b`",
        default=numpy.array([0.0]),
        range=basic.Range(lo=0.0, hi=30.0, step=0.1), #TODO need to define the range
        doc="""Inhibitory adaptation current increment [pA]""",
        order=11)

    a_i = arrays.FloatArray(
        label=":math:`Inhibitory a`",
        default=numpy.array([0.0]),
        range=basic.Range(lo=0.0, hi=20.0, step=0.1), #TODO need to define the range
        doc="""Excitatory adaptation conductance [nS]""",
        order=12)

    E_e = arrays.FloatArray(
        label=r":math:`E_e`",
        default=numpy.array([0.0]),
        range=basic.Range(lo=0.0, hi=60., step=0.01), #TODO need to define the range
        doc="""excitatory reversal potential [mV]""",
        order=13)

    E_i = arrays.FloatArray(
        label=":math:`E_i`",
        default=numpy.array([-80.0]),
        range=basic.Range(lo=-90.0, hi=-70.0, step=1.0), #TODO need to define the range
        doc="""inhibitory reversal potential [mV]""",
        order=14)

    Q_e = arrays.FloatArray(
        label=r":math:`Q_i`",
        default=numpy.array([1.0]),
        range=basic.Range(lo=0.0, hi=2.0, step=0.01), #TODO need to define the range
        doc="""excitatory quantal conductance [nS]""",
        order=15)

    Q_i = arrays.FloatArray(
        label=r":math:`Q_i`",
        default=numpy.array([5.0]),
        range=basic.Range(lo=0.0, hi=6.0, step=0.1), #TODO need to define the range
        doc="""inhibitory quantal conductance [nS]""",
        order=16)

    tau_e = arrays.FloatArray(
        label=":math:`\tau_e`",
        default=numpy.array([5.0]),
        range=basic.Range(lo=1.0, hi=6.0, step=1.0), #TODO need to define the range
        doc="""excitatory decay [ms]""",
        order=17)

    tau_i = arrays.FloatArray(
        label=":math:`\tau_i`",
        default=numpy.array([5.0]),
        range=basic.Range(lo=0.5, hi=2.0, step=0.01), #TODO need to define the range
        doc="""inhibitory decay [ms]""",
        order=18)

    N_tot = arrays.IntegerArray(
        label=":math:`N_{tot}`",
        default=numpy.array([10000]),
        range=basic.Range(lo=1000, hi=20000, step=1000), #TODO need to define the range
        doc="""cell number""",
        order=19)

    p_connect = arrays.FloatArray(
        label=":math:`\epsilon`",
        default=numpy.array([0.05]),
        range=basic.Range(lo=0.0, hi=1.0, step=0.01), #TODO need to define the range
        doc="""connectivity probability""",
        order=20)

    g = arrays.FloatArray(
        label=":math:`g`",
        default=numpy.array([0.2]),
        range=basic.Range(lo=0.0, hi=1.0, step=0.01), #TODO need to define the range
        doc="""fraction of inhibitory cells""",
        order=21)

    nu_ext = arrays.FloatArray(
        label=":math:`\nu_e^{drive}`",
        default=numpy.array([4.0]),
        range=basic.Range(lo=0.0, hi=20.0, step=0.01), #TODO need to define the range
        doc="""external drive""",
        order=22)

    T = arrays.FloatArray(
        label=":math:`T`",
        default=numpy.array([5.0]),
        range=basic.Range(lo=0.0, hi=20.0, step=0.1), #TODO need to define the range
        doc="""time scale of describing network activity""",
        order=22)

    P_e = arrays.IndexArray(
        label= ":math:`P_e`", #TODO need to check the size of the array when it's used
        required=False,
        file_storage=core.FILE_STORAGE_NONE, #TODO need to verified if the value is correct or not
        # default= numpy.array([]), #TODO need to define the fault value
        doc="""Polynome of excitatory phenomenological threshold (order 10)""",
        order=23)

    P_i = arrays.IndexArray(
        label=":math:`P_i`", #TODO need to check the size of the array when it's used
        required=False,
        file_storage=core.FILE_STORAGE_NONE, #TODO need to verified if the value is correct or not
        # default= numpy.array([]), #TODO need to define the fault value
        doc="""Polynome of inhibitory phenomenological threshold (order 10)""",
        order=24)


    # Used for phase-plane axis ranges and to bound random initial() conditions.
    state_variable_range = basic.Dict(
        label="State Variable ranges [lo, hi]",
        default={"E": numpy.array([0.0, 100.0]), #TODO need to define the range
                 "I": numpy.array([0.0, 100.0])}, #TODO need to define the range
        doc="""The values for each state-variable should be set to encompass
        the expected dynamic range of that state-variable for the current
        parameters, it is used as a mechanism for bounding random inital
        conditions when the simulation isn't started from an explicit history,
        it is also provides the default range of phase-plane plots.\n
        E: firing rate of excitatory population\n
        I: firing rate of inhibitory population
        """,
        order=23)

    variables_of_interest = basic.Enumerate(
        label="Variables watched by Monitors",
        options=["E", "I"], #TODO need to define the recordable variable
        default=["E"], #TODO need to define the default recordable variables
        select_multiple=True,
        doc="""This represents the default state-variables of this Model to be
               monitored. It can be overridden for each Monitor if desired. The
               corresponding state-variable indices for this model are :math:`E = 0`
               and :math:`I = 1`.""",
        order=24)

    state_variables = 'E I'.split()
    _nvar = 2
    cvar = numpy.array([0, 1], dtype=numpy.int32)

    def dfun(self, state_variables, coupling, local_coupling=0.0):
        r"""

        .. math::
            T \dot{\nu_\mu} &= -F_\mu(\nu_e,\nu_i) + \nu_\mu ,\all\mu\in\{e,i\}

        """
        #this variable can be internal varible for optimization
        TF = Transfer_Function(vars(self),self.P_e,self.P_i)

        E = state_variables[0, :]
        I = state_variables[1, :]
        derivative = numpy.empty_like(state_variables)

        # long-range coupling
        c_0 = coupling[0, :]

        # short-range (local) coupling
        lc_E = local_coupling * E
        lc_I = local_coupling * I


        #equation is inspired from github of Zerlaut :
        # https://bitbucket.org/yzerlaut/mean_field_for_multi_input_integration
        # ./mean_field/master_equation.py
        # Excitatory firing rate derivation
        derivative[0] = 1./self.T*(TF.excitatory(E+c_0+lc_E, I+lc_I)-E) #TODO : need to check the accuracy of the equation
        # Inhibitory firing rate derivation
        derivative[1] = 1./self.T*(TF.inhibitory(E+lc_E, I+lc_I)-I) #TODO : need to check the accuracy of the equation

        return derivative