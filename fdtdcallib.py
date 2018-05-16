import ctypes
from ctypes import c_int

import numpy as np
import numpy.ctypeslib as npct


# python 面向对象编程语言
class fdtdcallib(object):
    # 对类进行声明
    def __init__(self):
        # 初始化函数，引入动态链接库
        self.libc=npct.load_library("fdtdcal", ".")
        pass
    def initcal(self,ex1,ey1,ez1,bx1,by1,bz1,ex2,ey2,ez2,bx2,by2,bz2):
        # 类方法
        # 将计算空间电磁场的十二个分量传入C模块
        libcd=self.libc
        array_1d_double = npct.ndpointer(dtype=np.float64, ndim=3, flags='CONTIGUOUS') #声明指针类型级变量数据结构
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
        # 开始调用函数
        libcd.initecal(ex1,ey1,ez1,ff,dd)
        libcd.initbcal(bx1,by1,bz1,ff,dd)
        libcd.inite2cal(ex2,ey2,ez2,ff,dd)
        libcd.initb2cal(bx2,by2,bz2,ff,dd)
        return 'OK'



    def initc(self,ca,cb,cp,cq):
        # 对电磁介质参数进行初始化
        libcd=self.libc
        array_1d_double = npct.ndpointer(dtype=np.float64, ndim=3, flags='CONTIGUOUS')
        libcd.initc.restype=c_int
        libcd.initc.argtypes=[array_1d_double,array_1d_double,array_1d_double,array_1d_double]
        libcd.initc(ca,cb,cp,cq)
        return "OK"


    def initarr(self,n):
        # 对空间步长进行初始化
        libcd=self.libc
        arr=np.array(n,dtype='float64')
        libcd.intarry.restype=None
        libcd.intarry.argtypes=[ctypes.c_double]
        libcd.intarry(arr)
        return 'OK'
    def callib(self,n):
        libcd=self.libc
        nn=np.array(n,dtype=c_int)
        libcd.cal.restype=c_int
        libcd.cal.argtypes=[c_int]
        libcd.cal(nn)
        return 'OK'

    pass
