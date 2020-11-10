from tvb.simulator.noise import Additive,Noise,RandomStream
from tvb.datatypes import arrays
import numpy
from tvb.basic.traits import types_basic as basic

class  Ornstein_Ulhenbeck_process (Additive):
    tau_OU = basic.Float(
        label="time scale of decay",
        required=True,
        default=1.0,
        doc="""The noise time scale """)
    mu = arrays.FloatArray(
        label=":math:`mu`",
        required=True,
        default=numpy.array([1.0]),
        doc="""Mean of noise noise"""
    )
    weights=arrays.FloatArray(
        label=":math:`mu`",
        required=True,
        default=numpy.array([0.0]),
        doc="""Mean of noise noise"""
    )
    _noise= None

    def configure_white(self,dt,shape):
        """
        Run base classes configure to setup traited attributes, then ensure that
        the ``random_stream`` attribute is properly configured.

        """
        self.dt = dt
        self._noise = 0.0
        self.mu = numpy.reshape(self.mu,(7,1,1))

    def generate(self, shape, lo=-1.0, hi=1.0):
        self._noise = self._noise - self.dt/self.tau_OU*self._noise+numpy.sqrt(self.dt)*self.random_stream.normal(size=shape)
        noise = self.mu + self.nsig * self._noise
        return noise

    def gfun(self, state_variables):
        """
            Drop noise in order to avoid negative frequency

        """
        # drop value for negative noise
        return self.weights*1e-3


class Poisson_noise(Noise):
    nsig = arrays.FloatArray(
        label=r"rate",
        required=True,
        default=numpy.array([0.0]),
        doc="""rate of neurons in Hz"""
    )
    weights= arrays.FloatArray(
        label=r"unit firing rate",
        required=True,
        default=numpy.array([0.0]),
        doc="""TODO"""
    )

    noise_seed = basic.Integer(
        default=42,
        doc="A random seed used to initialise the random_stream if it is missing."
    )

    random_stream = RandomStream(
        required=False,
        label="Random Stream",
        doc="An instance of numpy's RandomState associated with this"
            "specific Noise object. Used when you need to resume a simulation from a state saved to disk"
    )

    def __init__(self, **kwargs):
        super(Poisson_noise, self).__init__(**kwargs)
        if self.random_stream is None:
            self.random_stream = numpy.random.RandomState(self.noise_seed)

        self.dt = None

    # For use if coloured
    def configure(self):
        """
        Run base classes configure to setup traited attributes, then ensure that
        the ``random_stream`` attribute is properly configured.

        """
        super(Poisson_noise, self).configure()
        # XXX: reseeding here will destroy a maybe carefully set random_stream!
        # self.random_stream.seed(self.noise_seed)

    def configure_white(self, dt, shape=None):
        """Set the time step (dt) of noise or integration time"""
        self.dt = dt

    def generate(self, shape, lo=0.0, hi=1.0):
        "Generate noise realization."
        lambda_poisson = self.nsig*self.dt*1e-3 # rate on the intervale dt (dt is in ms)
        noise = self.random_stream.poisson(lam=lambda_poisson,size=shape) # firing rate in Khz of the poisson generator
        # print(numpy.max(noise))
        noise[3]=numpy.sqrt(noise[0])
        return noise*self.weights

    def gfun(self, state_variables):
        r"""
        Linear additive noise, thus it ignores the state_variables.

        .. math::
            g(x) = \sqrt{2D}

        """
        g_x = numpy.ones(state_variables.shape)
        return g_x