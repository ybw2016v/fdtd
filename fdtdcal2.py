
# coding: utf-8

# In[2]:


import numpy as np
import numpy.ctypeslib as npct
import ctypes
from ctypes import c_int


# In[3]:


libcd = npct.load_library("fdtdcal", ".")


# In[4]:


def initcal(ex1,ey1,ez1,bx1,by1,bz1,ex2,ey2,ez2,bx2,by2,bz2,libcd):
    array_1d_double = npct.ndpointer(dtype=np.float64, ndim=3, flags='CONTIGUOUS')
    ff=np.array(ex1.strides,dtype='int32')
    dd=np.array(np.shape(ex1),dtype='int32')
    libcd.initecal.restype = c_int
    libcd.initbcal.restype = c_int
    libcd.inite2cal.restype = c_int
    libcd.initb2cal.restype = c_int
    libcd.initecal.argtypes=[array_1d_double,array_1d_double,array_1d_double,npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS'),npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS')]
    libcd.initbcal.argtypes=[array_1d_double,array_1d_double,array_1d_double,npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS'),npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS')]
    libcd.inite2cal.argtypes=[array_1d_double,array_1d_double,array_1d_double,npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS'),npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS')]
    libcd.initb2cal.argtypes=[array_1d_double,array_1d_double,array_1d_double,npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS'),npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS')]

    libcd.initecal(ex1,ey1,ez1,ff,dd)
    libcd.initbcal(bx1,by1,bz1,ff,dd)
    libcd.inite2cal(ex2,ey2,ez2,ff,dd)
    libcd.initb2cal(bx2,by2,bz2,ff,dd)
    return 'OK'


# In[5]:


def initarr(n,libcd):
    arr=np.array(n,dtype='float64')
    libcd.intarry.restype=None
    libcd.intarry.argtypes=[ctypes.c_double]
    libcd.intarry(arr)
    return 'OK'
    


# In[6]:


N=32
ex1=np.zeros([N,N,N])
ey1=np.zeros([N,N,N])
ez1=np.zeros([N,N,N])
ex2=np.zeros([N,N,N])
ey2=np.zeros([N,N,N])
ez2=np.zeros([N,N,N])
bx1=np.zeros([N,N,N])
by1=np.zeros([N,N,N])
bz1=np.zeros([N,N,N])
bx2=np.zeros([N,N,N])
by2=np.zeros([N,N,N])
bz2=np.zeros([N,N,N])
ex1[int(N/2),int(N/2),int(N/2)]=1


# In[7]:


def callib(n,libcd):
    nn=np.array(n,dtype=c_int)
    libcd.cal.restype=c_int
    libcd.cal.argtypes=[c_int]
    libcd.cal(nn)
    return 'OK'


# In[8]:


initarr(0.5,libcd)


# In[9]:


initcal(ex1,ey1,ez1,bx1,by1,bz1,ex2,ey2,ez2,bx2,by2,bz2,libcd)


# In[ ]:


callib(1,libcd)


# # In[ ]:


# ex.dtype





# # In[ ]:


# array_1d_double = npct.ndpointer(dtype=np.float64, ndim=3, flags='CONTIGUOUS')


# # In[ ]:


# ff=np.array(e1.strides,dtype='int32')


# # In[ ]:


# dd=np.array(np.shape(e1),dtype='int32')


# # In[ ]:


# libcd.initcal.restype = c_int


# In[ ]:


#libcd.initcal.argtypes=[array_1d_double,array_1d_double,array_1d_double,npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS'),npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS')]

