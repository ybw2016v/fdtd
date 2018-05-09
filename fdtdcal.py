
# coding: utf-8

# In[1]:


import numpy as np
import numpy.ctypeslib as npct
from ctypes import c_int


# In[2]:


ex=np.arange(1,101,1,dtype='float64')


# In[3]:


e1=ex.reshape(-1,2,2)


# In[4]:


ex.dtype


# In[5]:


e1


# In[6]:


libcd = npct.load_library("fdtdcal", ".")


# In[7]:


array_1d_double = npct.ndpointer(dtype=np.float, ndim=3, flags='CONTIGUOUS')


# In[8]:


ff=np.array(e1.strides,dtype='int32')


# In[9]:


dd=np.array(np.shape(e1),dtype='int32')


# In[10]:


libcd.initcal.restype = c_int


# In[11]:


libcd.initcal.argtypes=[array_1d_double,array_1d_double,array_1d_double,npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS'),npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS')]


# In[12]:


libcd.initcal(e1,e1,e1,ff,dd)

print(e1)