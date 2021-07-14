#!/usr/bin/env python2

import numpy as np
import matplotlib.pyplot as plt

omega0=10.0*2*np.pi
xi=0.1
M=1.
K=M*omega0**2
C=M*2*xi*omega0

omega = [2*2.0*np.pi]
while 1.06*omega[-1]<(25*2.0*np.pi):
    omega.append(1.06*omega[-1])
omega=np.array(omega)
X = np.abs(1.0/(-M*omega**2+C*1j*omega+K))

MAT1=np.genfromtxt("harmonicExcitationElem.usr")
omegambdyn1 = MAT1[:,2]
Xmbdyn1 = MAT1[:,4]

MAT3=np.genfromtxt("harmonicExcitationElem3.usr")
omegambdynMAT3 = MAT3[:,2]
omegambdyn3=np.unique(omegambdynMAT3)
Xmbdyn3 = 0*omegambdyn3
for i in range(len(omegambdyn3)):
    Xmbdyn3[i]=MAT3[np.argwhere(omegambdynMAT3==omegambdyn3[i])[-1],4]

MAT4=np.genfromtxt("harmonicExcitationElem4.usr")
omegambdynMAT4 = MAT4[:,2]
omegambdyn4=np.unique(omegambdynMAT4)
Xmbdyn4 = 0*omegambdyn4
for i in range(len(omegambdyn4)):
    Xmbdyn4[i]=1.0/MAT4[np.argwhere(omegambdynMAT4==omegambdyn4[i])[-1],-2]

MAT5=np.genfromtxt("harmonicExcitationElem5.usr")
omegambdynMAT5 = MAT5[:,2]
omegambdyn5=np.unique(omegambdynMAT5)
Xmbdyn5 = 0*omegambdyn5
for i in range(len(omegambdyn5)):
    Xmbdyn5[i]=MAT5[np.argwhere(omegambdynMAT5==omegambdyn5[i])[-1],4]

fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(omega/(2*np.pi),X,"-o",fillstyle='none',label="theoretical",)
ax.plot(omegambdyn1/(2*np.pi),Xmbdyn1,"-d",fillstyle='none',label="harmonicExcitationElem")
ax.plot(omegambdyn3/(2*np.pi),Xmbdyn3,"-s",fillstyle='none',label="harmonicExcitationElem3")
ax.plot(omegambdyn4/(2*np.pi),Xmbdyn4,"-x",fillstyle='none',label="harmonicExcitationElem4")
ax.plot(omegambdyn5/(2*np.pi),Xmbdyn5,"-*",fillstyle='none',label="harmonicExcitationElem5")
ax.set_xlabel("f, Hz")
ax.set_ylabel("transfer function, g/N")
ax.grid(True)
ax.set_yscale("log")
ax.legend(loc=0)
fig.savefig("harmonicExcitationElemResults.pdf",bbox_inches='tight')
