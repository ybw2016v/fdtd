#!/usr/bin/env python3
import math

import numpy as np

from fdtdcallib import *


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
N=128
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
ca=np.ones([N,N,N])*1 #(((1-(thegam*dt)/(2*eps))/(1+(thegam*dt)/(2*eps))))
cb=np.ones([N,N,N])*(dt/(eps*d)) #((dt/(d*eps))/(1+thegam*dt/(2*eps)))
cp=np.ones([N,N,N])*1 #(((1-(miy*dt)/(2*miy*miu))/(1+(miy*dt)/(2*miy*miu))))
cq=np.ones([N,N,N])*(dt/(miu*d)) #((dt/(d*miu))/(1+miy*dt/(2*miu)))

fdtd=fdtdcallib()

fdtd.initcal(ex1,ey1,ez1,bx1,by1,bz1,ex2,ey2,ez2,bx2,by2,bz2)

fdtd.initc(ca,cb,cp,cq)

for iii in range(0,20):
    ez1[32,32,12:24]=1*math.sin(2*omega*dt*iii)
#     ey1[int(N/2),int(N/2),int(N/2)]=0
#     ez1[int(N/2),int(N/2),int(N/2)]=0
    print('*')
    print(ez1[32,32,12:24])
    fdtd.callib(1)
    print(ez1[32,32,12:24])
    print('\n')
np.save('ex',ex1)
