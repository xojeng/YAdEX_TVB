from tvb.simulator.noise import Noise
from tvb.datatypes import arrays
import numpy
from tvb.basic.traits import types_basic as basic

class  Ornstein_Ulhenbeck_process (Noise):

    nsig = arrays.FloatArray(
        configurable_noise=True,
        label=":math:`D`",
        required=True,
        default=numpy.array([1.0]), range=basic.Range(lo=0.0, hi=10.0, step=0.1),
        order=1,
        doc=""" the noise amplitude \n
        NOTE:
        Sensible values are typically ~<< 1% of the dynamic range of a Model's
        state variables.""")
    T = arrays.FloatArray(
        configurable_noise=True,
        label=":math:`D`",
        required=True,
        default=numpy.array([1.0]), range=basic.Range(lo=0.0, hi=10.0, step=0.1),
        order=1,
        doc="""fast time scale of the process """)
    _noise =0.0
        
    def generate(self, shape, lo=-1.0, hi=1.0):
        """
            Generate a noise for all the shape : Ornstein-Uhlenbeck process
        """
        derive = self._noise/self.T*self.dt + numpy.sqrt(self.dt) * self.random_stream.normal(size=shape)
        self._noise += self.dt*self.nsig*derive
        return self._noise

    def gfun(self, state_variables):
        """
            Drop noise in order to avoid negative frequency

        """
        g_x = numpy.ones(state_variables.shape)
        check = self._noise + state_variables
        # drop value for negative noise
        g_x[:,numpy.where(check[0,:,:] <= 0.0)] = 0.0
        g_x[:,numpy.where(check[1,:,:] <= 0.0)] = 0.0
        return g_x
