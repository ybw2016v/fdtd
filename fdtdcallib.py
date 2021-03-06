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
    def initcal(self,ex1,ey1,ez1,bx1,by1,bz1):
        # 类方法
        # 将计算空间电磁场的十二个分量传入C模块
        libcd=self.libc
        array_1d_double = npct.ndpointer(dtype=np.float64, ndim=3, flags='CONTIGUOUS') #声明指针类型级变量数据结构
        ff=np.array(ex1.strides,dtype='int32')
        dd=np.array(np.shape(ex1),dtype='int32')
        libcd.initecal.restype = c_int
        libcd.initbcal.restype = c_int
        # libcd.inite2cal.restype = c_int
        # libcd.initb2cal.restype = c_int
        libcd.initecal.argtypes=[array_1d_double,array_1d_double,array_1d_double,npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS'),npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS')]
        libcd.initbcal.argtypes=[array_1d_double,array_1d_double,array_1d_double,npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS'),npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS')]
        # libcd.inite2cal.argtypes=[array_1d_double,array_1d_double,array_1d_double,npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS'),npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS')]
        # libcd.initb2cal.argtypes=[array_1d_double,array_1d_double,array_1d_double,npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS'),npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS')]
        # 开始调用函数
        libcd.initecal(ex1,ey1,ez1,ff,dd)
        libcd.initbcal(bx1,by1,bz1,ff,dd)
        # libcd.inite2cal(ex2,ey2,ez2,ff,dd)
        # libcd.initb2cal(bx2,by2,bz2,ff,dd)
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
    def initmur(self,c,dt,dx):
        libcd=self.libc
        cc=np.array(c,dtype='float64')
        ddt=np.array(dt,dtype='float64')
        ddx=np.array(dx,dtype='float64')
        libcd.initmur1.restype=c_int
        libcd.initmur1.argtypes=[ctypes.c_double,ctypes.c_double,ctypes.c_double]
        libcd.initmur1(cc,ddt,ddx)
        return 'OK'
    def calmur1lib(self,n):
        libcd=self.libc
        nn=np.array(n,dtype=c_int)
        libcd.calmur1.restype=c_int
        libcd.calmur1.argtypes=[c_int]
        libcd.calmur1(nn)
        return 'OK'
    pass
    def calmur2lib(self,n):
        libcd=self.libc
        nn=np.array(n,dtype=c_int)
        libcd.calmur2.restype=c_int
        libcd.calmur2.argtypes=[c_int]
        libcd.calmur2(nn)
        return 'OK'
    pass
    def initmur2(self,ex1,ey1,ez1):
        libcd=self.libc
        array_1d_double = npct.ndpointer(dtype=np.float64, ndim=3, flags='CONTIGUOUS') #声明指针类型级变量数据结构
        ff=np.array(ex1.strides,dtype='int32')
        dd=np.array(np.shape(ex1),dtype='int32')
        libcd.initmur2.argtypes=[array_1d_double,array_1d_double,array_1d_double,npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS'),npct.ndpointer(dtype=c_int,ndim=1,flags='CONTIGUOUS')]
        libcd.initmur2.restype = c_int
        libcd.initmur2(ex1,ey1,ez1,ff,dd)
        return 'OK'
    pass
