# from progressbar import * 
import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import pyplot as plt
e=np.load('ex.npy')
plt.figure(1)
plt.imshow(e[:,:,128])
plt.colorbar()
plt.savefig('dog.jpg')
plt.close(1)
plt.figure(2)
k=np.load('k.npy')
plt.plot(k)
plt.savefig('kk.jpg')
plt.close(2)
# pbar = ProgressBar()
for iii in range(1,26):
    plt.figure(iii)
    adog=np.load('r'+str(iii)+'.npy')
    plt.imshow(adog[100:140,90:150,105],vmin=-1, vmax=1)
    plt.colorbar()
    plt.savefig('r'+str(iii)+'.jpg')
    plt.close(iii)
    pass
