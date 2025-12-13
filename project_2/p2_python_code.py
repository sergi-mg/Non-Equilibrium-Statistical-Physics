# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 09:15:15 2025

@author: Sergi Martinez Galindo
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import random as rand
from tqdm import tqdm
#%%
def calc_reg(x,y):
    import numpy as np
    #y=mx+b
    N=np.size(x)
    #calcul mitjanes
    xmed=np.sum(x)/N
    x2med=np.sum(x**2)/N
    ymed=np.sum(y)/N
    y2med=np.sum(y**2)/N
    xymed=np.sum(x*y)/N
    #calcul sigmes
    sigx2=x2med-xmed**2
    sigy2=y2med-ymed**2
    sigxy=xymed-xmed*ymed
    #calcul regressio
    m=sigxy/sigx2
    b=ymed-m*xmed
    r=sigxy/((sigx2*sigy2)**0.5)
    #calcul incerteses
    dyreg=((sigy2*(1-r**2)*N)/(N-2))**0.5
    dm=dyreg/(N*sigx2)**0.5
    db=dyreg*(x2med/(N*sigx2))**0.5
    #retornem els valors
    return m,dm,b,db,r

def calc_reg_no_ord(x,y):
    import numpy as np
    #y=kx
    N=np.size(x)
    sumx2=np.sum(x**2)
    sumxy=np.sum(x*y)
    k=sumxy/sumx2
    sumy_kx=np.sum((y-k*x)**2)
    dyreg=((sumy_kx)/(N-1))**0.5
    dk=dyreg/(sumx2)**0.5
    return k, dk

def xifres_escalar(valor,error,exp_max,exp_min):
    valor_c=0
    error_c=0
    for k in range(exp_max,exp_min,-1):
        if error>=2*10**k:
            valor_c=round(valor,-k)
            error_c=round(error,-k)
            break
    return valor_c,error_c
#%%
#Section 2
dades_1=np.loadtxt("delta_t.dat")
#%%
dades_1_tau1=dades_1[:2001]
dades_1_tau2=dades_1[2001:4002]
dades_1_tau3=dades_1[4002:]
dades_list=[dades_1_tau1,dades_1_tau2,dades_1_tau3]
tau_legend=[1,10,100]
D=15
#plot
cmap = plt.cm.viridis
plt.figure()
col=-0.4
for i in range(3):
    print("tau=",tau_legend[i])
    col+=0.4
    result=dades_list[i]
    #plot
    x=result[1:,0]
    f=2*D*(x+tau_legend[i]*(np.exp(-x/tau_legend[i])-1))
    
    plt.plot(result[:,0], result[:,1],linestyle="none",\
                 color=cmap(col),marker="o",markersize=3,\
                     label=r"$\tau="+str(tau_legend[i])+"$")
    plt.plot(x,f,linestyle="dashed",color=cmap(col+0.15),\
             label=r"Theoretical: $\tau="+str(tau_legend[i])+"$")

x=np.arange(0.05,10000,0.05)
f=2*D*x
plt.plot(x,f,linestyle="dashed",color="black",\
         label=r"Theoretical: Passive Brownian")        
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$t$",fontsize=12)
plt.ylabel("$\Delta^2$",fontsize=12)
plt.xlim((0.05,10000))
plt.legend(loc="upper center",bbox_to_anchor=(0.5, 1.32), ncol=2 ,fontsize=12)
plt.grid(which="major", linestyle='-', linewidth=0.7)
plt.grid(which="minor", linestyle='--', linewidth=0.3)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.savefig("Im2.pdf", bbox_inches='tight')
plt.show()

#%%
#Section 3
dades=np.loadtxt("correlacio.dat")
#%%
dades_tau1=dades[:2003001]
dades_tau2=dades[2003001:4006002]
dades_tau3=dades[4006002:]
datasets=[dades_tau1,dades_tau2,dades_tau3]
cmap = plt.cm.viridis
col=-0.4
tau_legend=[1,10,100]
plt.figure()
for i in range(3):
    print("tau=",tau_legend[i])
    col+=0.4
    delta_t = datasets[i][:, 0]
    corr = datasets[i][:, 1] #correlation function
    
    # delta t values and inverse index
    vals, inv = np.unique(delta_t, return_inverse=True)
    
    # pmean value for each delta t
    corr_mean = np.bincount(inv, weights=corr) / np.bincount(inv)
    
    # final result
    result = np.column_stack((vals, corr_mean))
    
    #linear regression
    a,da,b,db,r=calc_reg(result[:81,0],np.log(result[:81,1]))
    #tau
    t,dt=xifres_escalar(-1/a,da/a**2,10,-10)
    print("Tau=",t,"+-",dt)
    #D
    D,dD=xifres_escalar(t*np.exp(b),np.exp(b)*(dt+t*db),10,-10)
    print("D=",D,"+-",dD)
    
    #regression in the graphic
    C,dC=xifres_escalar(np.exp(b),np.exp(b)*(db),10,-10)
    t_reg=result[:,0]
    
    #plot
    plt.plot(result[:,0], result[:,1]/C,linestyle="none",\
                 color=cmap(col),marker="o",markersize=3,\
                     label=r"$\tau="+str(tau_legend[i])+"$")
    plt.plot(t_reg,np.exp(-t_reg/t),linestyle="dashed",color=cmap(col+0.15),\
             label=r"Linear regression: $\tau="+str(tau_legend[i])+"$")
        
#plt.yscale("log")
plt.xlim((-10,1000))
plt.xlabel(r"$|t-t'|$",fontsize=12)
plt.ylabel(r"$\langle\xi(t)\xi(t')\rangle\cdot\tau/D$",fontsize=12)
plt.legend(loc="upper center",bbox_to_anchor=(0.5, 1.3), ncol=2 ,fontsize=12)
plt.grid(which="major", linestyle='-', linewidth=0.7)
plt.grid(which="minor", linestyle='--', linewidth=0.3)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.savefig("Im3.pdf", bbox_inches='tight')

plt.show()

#%%
dades_2=np.loadtxt("d_force_t.dat")
#%%
f_data=[]
for i in range(5):
    f_data.append(dades_2[2001*i:2001*(i+1)])
    
#plot
fig1, ax1 = plt.subplots() #MSD
fig2, ax2 = plt.subplots() #displacement
f_legend=[0,0.001,0.01,0.1,1]
cmap = plt.cm.viridis
col=-0.2
for i in range(4,-1,-1):
    print("f=",f_legend[i])
    col+=0.2
    result=f_data[i]
    
    ax1.plot(result[:,0], result[:,1],linestyle="none",\
                 color=cmap(col),marker="o",markersize=1.5,\
                     label=r"$f="+str(f_legend[i])+"$")
    ax2.plot(result[:,0], result[:,2],linestyle="none",\
                 color=cmap(col),marker="o",markersize=1.5,\
                     label=r"$f="+str(f_legend[i])+"$")
        
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel(r"$t$",fontsize=12)
ax1.set_ylabel(r"$\Delta^2$",fontsize=12)
ax1.legend(loc="upper center",bbox_to_anchor=(0.5, 1.25), ncol=2 ,fontsize=12)
ax1.grid(which="major", linestyle='-', linewidth=0.7)
ax1.grid(which="minor", linestyle='--', linewidth=0.3)
ax1.tick_params(axis='both', which='major', labelsize=12)
fig1.savefig("Im4_1.pdf", bbox_inches='tight')
fig1.show()
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel(r"$t$",fontsize=12)
ax2.set_ylabel(r"$d$",fontsize=12)
ax2.legend(loc="upper center",bbox_to_anchor=(0.5, 1.25), ncol=2 ,fontsize=12)
ax2.grid(which="major", linestyle='-', linewidth=0.7)
ax2.grid(which="minor", linestyle='--', linewidth=0.3)
ax2.tick_params(axis='both', which='major', labelsize=12)
fig2.savefig("Im4_2.pdf", bbox_inches='tight')
fig2.show()

#%%
dades_3=np.loadtxt("d_force_t_2.dat")
f_t_data=[]
for i in range(30):
    f_t_data.append(dades_3[2001*i:2001*(i+1)])
    
f_val=np.array([0,0.001,0.01,0.1,1])
tau_val=np.array([1,20,40,60,80,100])

D=np.zeros((6,5,2))
mu=np.zeros((6,4,2))

data_index=-1
for i in range(len(tau_val)):
    tau=tau_val[i]
    for j in range(len(f_val)):
        data_index+=1
        f=f_val[j]
        data=f_t_data[data_index]
        
        y_MSD=data[:-250,1]
        y_d=data[:-250,2]
        x_t=data[:-250,0]
        
        #linear regression
        a1,da1=calc_reg_no_ord(x_t,y_MSD)
        D_ij,dD_ij=xifres_escalar(a1/2,da1/2,5,-7)
        D[i,j,:]=np.array([D_ij,dD_ij])
        
        if f>0:
            a2,da2=calc_reg_no_ord(x_t,y_d)
            mu_ij,dmu_ij=xifres_escalar(a2/f,da2/f,5,-7)
            mu[i,j-1,:]=np.array([mu_ij,dmu_ij])
        

fig1, ax1 = plt.subplots() #D(tau)
fig2, ax2 = plt.subplots() #mu(tau)       
cmap = plt.cm.viridis
col=-0.2
for i in range(len(f_val)):
    print("f=",f_val[i])
    f=f_val[i]
    col+=0.2
    result=f_data[i]
    
    ax1.errorbar(tau_val, D[:,i,0],xerr=0,yerr=D[:,i,1],linestyle="dashed",\
                 color=cmap(col),marker="o",markersize=3,\
                    label=r"$f="+str(f_legend[i])+"$")
    if f>0:
        ax2.errorbar(tau_val, mu[:,i-1,0],xerr=0,yerr=mu[:,i-1,1],linestyle="dashed",\
                     color=cmap(col),marker="o",markersize=3,\
                        label=r"$f="+str(f_legend[i])+"$")
        
#ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel(r"$\tau$",fontsize=12)
ax1.set_ylabel(r"$D$",fontsize=12)
ax1.legend(loc="upper center",bbox_to_anchor=(0.5, 1.25), ncol=2 ,fontsize=12)
ax1.grid(which="major", linestyle='-', linewidth=0.7)
ax1.grid(which="minor", linestyle='--', linewidth=0.3)
ax1.tick_params(axis='both', which='major', labelsize=12)
fig1.savefig("Im5_1.pdf", bbox_inches='tight')
fig1.show()
#ax2.set_xscale("log")
#ax2.set_yscale("log")
ax2.set_xlabel(r"$\tau$",fontsize=12)
ax2.set_ylabel(r"$\mu$",fontsize=12)
ax2.legend(loc="upper center",bbox_to_anchor=(0.5, 1.2), ncol=2 ,fontsize=12)
ax2.grid(which="major", linestyle='-', linewidth=0.7)
ax2.grid(which="minor", linestyle='--', linewidth=0.3)
ax2.tick_params(axis='both', which='major', labelsize=12)
fig2.savefig("Im5_2.pdf", bbox_inches='tight')
fig2.show()
