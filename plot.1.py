# from progressbar import * 
import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import pyplot as plt
e=np.load('ex.npy')
plt.figure(1)
plt.imshow(e[:,:,16])
plt.colorbar()
plt.savefig('dog.jpg')
plt.close(1)
# plt.figure(2)
# # k=np.load('k.npy')
# plt.plot(k)
# plt.savefig('kk.jpg')
# plt.close(2)
# pbar = ProgressBar()
for iii in range(1,50):
    plt.figure(iii)
    adog=np.load('s'+str(iii)+'.npy')
    plt.imshow(adog[:,16,:])
    plt.colorbar()
    plt.savefig('s'+str(iii)+'.jpg')
    plt.close(iii)
    pass
