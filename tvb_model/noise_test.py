import noise
import numpy as np

np.random.seed(10)
noise_test = noise.Ornstein_Ulhenbeck_process(nsig=np.array([[1.0e-6],[1.0e-6],[0.0],[0.0],[0.0],[0.0],[0.0]]))
noise_test.dt = 0.5
nb_test = 800000*2
nb_region = 76
value=np.zeros((nb_test,7,nb_region))
for i in range(nb_test):
    value[i,:,:] = np.squeeze(noise_test.generate((7,nb_region,1))*noise_test.gfun(np.ones((7,nb_region,1))))
    print('min noise :'+str(np.min(value[i,:,:])))
    print('max noise :'+str(np.max(value[i,:,:])))

np.save('./value.npy',value)