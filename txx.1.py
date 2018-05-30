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
N=256
ex1=np.zeros([N,N,N])
ey1=np.zeros([N,N,N])
ez1=np.zeros([N,N,N])
ex2=np.zeros([N,N,N])
ey2=np.zeros([N,N,N])
ez2=np.zeros([N,N,N])
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
ca[80,128,90:150]=1.5*0.2
cb[80,128,90:150]=217.5053478890454*0.005
fdtd=fdtdcallib()

fdtd.initcal(ex1,ey1,ez1,bx1,by1,bz1)

fdtd.initc(ca,cb,cp,cq)
pbar = ProgressBar()

num=3500

fdtd.initmur(299792458.0,2.309401076758503e-12,0.0011991698319999999)#(1,0.5,1)#
fdtd.initmur2(ex2,ey2,ez2)
k=np.zeros([num,5])
for iii in pbar(range(0,num)):
    
#     ey1[int(N/2),int(N/2),int(N/2)]=0
#     ez1[int(N/2),int(N/2),int(N/2)]=0
    # print(iii)
    # print(ez1[16,16,11:24])
    # pbar.update(int((iii / (500 - 1)) * 100))
    # fdtd.callib(1)
    fdtd.calmur2lib(1)
    ez1[128,128,100:140]=1*math.sin(1*0.01451039491387374*iii)#2*omega*dt*
    # if iii<=10000:
    #     ez1[32,32,30:34]=1*math.sin(5*0.01451039491387374*iii)#2*omega*dt*
    # else:
    #     ez1[32,32,30:34]=0

    if iii%100==0:
        np.save('s'+str(int(iii/100)),ez1)
    # print(ez1[.016,16,11:24])
    # print('\n')
    k[iii,0]=ez1[160,20,128]
    k[iii,1]=ez1[160,128,200]
    k[iii,2]=ez1[160,180,50]
    k[iii,3]=ez1[160,200,50]
    k[iii,4]=ez1[160,20,50]
# pbar.finish()
E=ex1**2+ey1**2+ez1**2
np.save('sex',E)
np.save('sk',k)
  
