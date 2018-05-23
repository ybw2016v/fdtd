import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import pyplot as plt
e=np.load('ex.npy')
plt.figure(1)
plt.imshow(e[:,:,32])
plt.colorbar()
plt.savefig('dog.jpg')

plt.figure(2)
k=np.load('k.npy')
plt.plot(k)
plt.savefig('kk.jpg')
