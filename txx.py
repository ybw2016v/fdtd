#!/usr/bin/env python3
import math

import numpy as np

from fdtdcallib import *

from progressbar import * 

miu=4*math.pi*10e-8
# 真空磁导率
eps=8.85*10e-13
# 绝对介电常数
thegam=0
# 电导率
miy=1
# 导磁率
d=3e-5
# 空间步长m
dt=1/((3e9)*math.sqrt(3)*d)
# 时间步长
FREQUENCY=10.0e9
# 信号频率
omega = 2.0*math.pi*FREQUENCY
N=64
ex1=np.zeros([N,N,N])
ey1=np.zeros([N,N,N])
ez1=np.zeros([N,N,N])
# ex2=np.zeros([N,N,N])
# ey2=np.zeros([N,N,N])
# ez2=np.zeros([N,N,N])
bx1=np.zeros([N,N,N])
by1=np.zeros([N,N,N])
bz1=np.zeros([N,N,N])
# bx2=np.zeros([N,N,N])
# by2=np.zeros([N,N,N])
# bz2=np.zeros([N,N,N])
ca=np.ones([N,N,N])*1 #(((1-(thegam*dt)/(2*eps))/(1+(thegam*dt)/(2*eps))))
cb=np.ones([N,N,N])*217.5053478890454*0.1 #*(dt/(eps*d))#((dt/(d*eps))/(1+thegam*dt/(2*eps)))
cp=np.ones([N,N,N])*1 #(((1-(miy*dt)/(2*miy*miu))/(1+(miy*dt)/(2*miy*miu))))
cq=np.ones([N,N,N])*0.0015325293679830556*0.1# *(dt/(miu*d))#((dt/(d*miu))/(1+miy*dt/(2*miu)))

fdtd=fdtdcallib()

fdtd.initcal(ex1,ey1,ez1,bx1,by1,bz1)

fdtd.initc(ca,cb,cp,cq)
pbar = ProgressBar()

num=500

fdtd.initmur(299792458.0,2.309401076758503e-12,0.0011991698319999999)

k=np.zeros(num)
for iii in pbar(range(0,num)):
    if iii<=100:
        ez1[10,32,31:33]=1*math.sin(1*0.01451039491387374*iii)#2*omega*dt*
    else:
        ez1[10,32,31:33]=0
#     ey1[int(N/2),int(N/2),int(N/2)]=0
#     ez1[int(N/2),int(N/2),int(N/2)]=0
    # print(iii)
    # print(ez1[16,16,11:24])
    # pbar.update(int((iii / (500 - 1)) * 100))
    # fdtd.callib(1)
    fdtd.calmur1lib(1)
    # print(ez1[16,16,11:24])
    # print('\n')
    k[iii]=ez1[16,17,15]
# pbar.finish()
E=ex1**2+ey1**2+ez1**2
np.save('ex',E)
np.save('k',k)
  
