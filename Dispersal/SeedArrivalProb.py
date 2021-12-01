#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 19:02:07 2021

@author: yanlanliu
"""

import numpy as np
import pandas as pd
import os
from scipy import ndimage
import mahotas


SHRUB_TYPE = 3
datapath ='../Input/'
datapath = '/global/scratch/yanlanliu/NGEE_data/Input/cleaned_code/'
lc0name = 'lc0'
gname, Bh, Bv = ('Bh04v01', 4, 1)
N, Nclm, dx = (6000,50,30)
SCALE_CLM = int(N/Nclm)


#======= Dispersal kernel parameters ========
theta_l = np.array([3900,0.5])
theta_s = np.array([600,1.5])
#============================================


#%%
def readLC1(gname,lc0name,fill_NA = 255):
    fname = datapath+gname+'.csv'
    if os.path.isfile(fname):
        df = pd.read_csv(fname)[['r','c',lc0name]]
        LC = df2img(df,lc0name,N,N,'r','c')
    else:
        LC = np.zeros([N,N],dtype=int)+fill_NA
    return LC

def f_add_padding(LC,Bh,Bv,lc0name,rmax):
    left = 'Bh'+str(Bh-1).zfill(2)+'v'+str(Bv).zfill(2)
    right = 'Bh'+str(Bh+1).zfill(2)+'v'+str(Bv).zfill(2)
    up = 'Bh'+str(Bh).zfill(2)+'v'+str(Bv-1).zfill(2)
    down = 'Bh'+str(Bh).zfill(2)+'v'+str(Bv+1).zfill(2)
    nw = 'Bh'+str(Bh-1).zfill(2)+'v'+str(Bv-1).zfill(2)
    ne = 'Bh'+str(Bh+1).zfill(2)+'v'+str(Bv-1).zfill(2)
    sw = 'Bh'+str(Bh-1).zfill(2)+'v'+str(Bv+1).zfill(2)
    se = 'Bh'+str(Bh+1).zfill(2)+'v'+str(Bv+1).zfill(2)
    row1 = np.concatenate([readLC1(nw,lc0name)[-rmax:,-rmax:],readLC1(up,lc0name)[-rmax:,:],readLC1(ne,lc0name)[-rmax:,:rmax]],axis=1)
    row2 = np.concatenate([readLC1(left,lc0name)[:,-rmax:],LC,readLC1(right,lc0name)[:,:rmax]],axis=1)
    row3 = np.concatenate([readLC1(sw,lc0name)[:rmax,-rmax:],readLC1(down,lc0name)[:rmax,:],readLC1(se,lc0name)[:rmax,:rmax]],axis=1)
    return np.concatenate([row1,row2,row3],axis=0)

def f_kernel(x,theta):
    return np.exp(-(x/theta[0])**theta[1])
    
def f_pr_seed(shb_presence,f_kernel,dist,theta):
    k = f_kernel(dist,theta)*dx # weight matrix
    Pr = ndimage.convolve(shb_presence,k,mode='constant',cval=0)
    rmax = int(dist.shape[0]/2)
    return Pr[rmax:-rmax,rmax:-rmax]/np.sum(k)


def caldist(rmax,trd=np.nan,theta=np.nan):
    dist = np.ones([rmax*2,rmax*2]); dist[rmax,rmax] = 0
    dist = np.sqrt(mahotas.distance(dist))*dx
    dist[dist==0] = 1 # to prevent nan in log-normal kernel
    if trd>0:
        disp_tmp = f_kernel(dist,theta)
        disp_tmp[disp_tmp<trd] = np.nan
        rmin,cmin = [min(tt) for tt in np.where(~np.isnan(disp_tmp))]
        rmax,cmax = [max(tt) for tt in np.where(~np.isnan(disp_tmp))]
        dist = dist[rmin:rmax,cmin:cmax]
    return dist

def df2img(df_lcc,vname,nr,nc,rname,cname):    
    r = np.arange(nr); c = np.arange(nc)
    ID_C, ID_R  =np.meshgrid(c,r)
    df0 = pd.DataFrame({rname:ID_R.flatten(), cname:ID_C.flatten()},index=np.arange(len(ID_R.flatten())))
    df_lcc = pd.merge(df0, df_lcc, how='left', on=[rname, cname])
    return df_lcc.pivot(index=rname, columns=cname, values=vname).values

#rmax_s, rmax_l = (2000, 100000) # maximum dispersal distance in m
rmax_s, rmax_l = (100,2000) # for testing (less computing time and memory)

rmax_s = int(rmax_s/dx); rmax_l = int(rmax_l/dx) # m to # of pixels

dist_s = caldist(rmax_s,trd=1e-3,theta=theta_s) #;plt.imshow(dist_s)
dist_l = caldist(rmax_l)
dist_zoom = ndimage.zoom(dist_l,1/24,order=1)


def f_dispersal(lc,lc0name,theta_s,theta_l):
    LCimg = df2img(lc,lc0name,N,N,'r','c')
    LCimg = f_add_padding(LCimg,Bh,Bv,lc0name,rmax_l) # useded when adjacent Bgrids are included
    LCimg = (LCimg==SHRUB_TYPE)*1
    # calculate long-distance dispersal
    Lzoom = ndimage.zoom(LCimg,1/24,order=1)
    Prl = f_pr_seed(Lzoom,f_kernel,dist_zoom,theta_l)
    Prl = ndimage.zoom(Prl,24,order=1)[lc['r'],lc['c']]
    # calculate short-distance dispersal
    LCimg = LCimg[rmax_l-rmax_s:rmax_s-rmax_l,rmax_l-rmax_s:rmax_s-rmax_l] 
    Prs = f_pr_seed(LCimg,f_kernel,dist_s,theta_s)
    Prs = Prs[lc['r'],lc['c']]
    Prl[np.isnan(Prl)] = 0; Prs[np.isnan(Prs)] = 0
    return Prl,Prs
    
#%%

fname = datapath+'LC_TOPO_'+gname+'.csv'
lc = pd.read_csv(fname)[['r','c',lc0name]]
lc['displ'],lc['disps'] = f_dispersal(lc,lc0name,theta_s,theta_l)

print(lc.head()) 
