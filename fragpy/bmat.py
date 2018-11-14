# Copyright 2018 The Fragpy Developers. All Rights Reserved.
# 
# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
# 
#   http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.
from fragpy.opt_param import param
from numpy import *
import numpy as np 
import numpy.linalg as linalg
from numpy.linalg import norm

'''
Forms Wilson's B Matrix 
'''

def BmatB(Q,int_coords): 

    intcoord = int_coords.intcoord 
    tmpclist = int_coords.tmpclist 
    anglist  = int_coords.anglist  
    dih      = int_coords.dih      
    hbonlist = int_coords.hbonlist 
    Hangle   = int_coords.Hangle   
    Hdih     = int_coords.Hdih     

    #---- Wilson's B matrix and its pseudo inverse ----
    bmatINT= np.zeros((intcoord,3*param.natom))
    bcount=0
    for i in tmpclist:
        dQbdxyz=bonval(Q[int(i[0])-1].q,Q[int(i[1])-1].q)
        for j in range(3):
            bmatINT[bcount][3*(int(i[0])-1)+j]=dQbdxyz[0][j]
            bmatINT[bcount][3*(int(i[1])-1)+j]=dQbdxyz[1][j]
        bcount+=1
        
    for i in anglist:
        dQadxyz=angval(Q[int(i[0])-1].q,Q[int(i[1])-1].q,Q[int(i[2])-1].q)
        for j in range(3):
            bmatINT[bcount][3*(int(i[0])-1)+j]=dQadxyz[0][j]
            bmatINT[bcount][3*(int(i[1])-1)+j]=dQadxyz[1][j]
            bmatINT[bcount][3*(int(i[2])-1)+j]=dQadxyz[2][j]
        bcount+=1

    for i in dih:
        dQddxyz=dihval(Q[int(i[0])-1].q,Q[int(i[1])-1].q,Q[int(i[2])-1].q,\
                       Q[int(i[3])-1].q)
        for j in range(3):
            bmatINT[bcount][3*(int(i[0])-1)+j]=dQddxyz[0][j]
            bmatINT[bcount][3*(int(i[1])-1)+j]=dQddxyz[1][j]
            bmatINT[bcount][3*(int(i[2])-1)+j]=dQddxyz[2][j]
            bmatINT[bcount][3*(int(i[3])-1)+j]=dQddxyz[3][j]
        bcount+=1
    
    for i in hbonlist:
        dQbhdxyz=bonval(Q[int(i[0])-1].q,Q[int(i[1])-1].q)
        for j in range(3):
            bmatINT[bcount][3*(int(i[0])-1)+j]=dQbhdxyz[0][j]
            bmatINT[bcount][3*(int(i[1])-1)+j]=dQbhdxyz[1][j]
        bcount+=1

    for i in Hangle:
        dQahdxyz=angval(Q[int(i[0])-1].q,Q[int(i[1])-1].q,Q[int(i[2])-1].q)
        for j in range(3):
            bmatINT[bcount][3*(int(i[0])-1)+j]=dQahdxyz[0][j]
            bmatINT[bcount][3*(int(i[1])-1)+j]=dQahdxyz[1][j]
            bmatINT[bcount][3*(int(i[2])-1)+j]=dQahdxyz[2][j]
        bcount+=1

    for i in Hdih:
        dQdhdxyz=dihval(Q[int(i[0])-1].q,Q[int(i[1])-1].q,Q[int(i[2])-1].q,\
                        Q[int(i[3])-1].q)
        for j in range(3):
            bmatINT[bcount][3*(int(i[0])-1)+j]=dQdhdxyz[0][j]
            bmatINT[bcount][3*(int(i[1])-1)+j]=dQdhdxyz[1][j]
            bmatINT[bcount][3*(int(i[2])-1)+j]=dQdhdxyz[2][j]
            bmatINT[bcount][3*(int(i[3])-1)+j]=dQdhdxyz[3][j]
        bcount+=1
    bcount = 0
    #-- pseudo inverse of B matrix --
    BinvINT=np.linalg.pinv(bmatINT,rcond=1e-15)

    return(bmatINT,BinvINT)

def PmatP(mat1,intcoord):
    G_mat = np.dot(np.array(mat1),np.array(mat1).T)
    G_inv = symm_matrix_inv(G_mat,intcoord)
    p_mat = np.dot(G_mat,G_inv)
    return(p_mat,G_inv)

def symm_matrix_inv(Alist1,dim):
    det = 1.0
    evals,A_evects = np.linalg.eigh(Alist1)

    for i in range(dim):
        det *= evals[i]

    A_inv = np.zeros((dim,dim))

    for i in range(dim):
        if fabs(evals[i]) > 1.0e-10:
            A_inv[i][i] = 1.0/evals[i]

    A_temp = np.dot(A_inv,A_evects.T)
    A_inv = np.dot(A_evects,A_temp)
    return A_inv

# ---> Bond
def bonval(m,n):
    u = subtract(m,n).tolist()
    r = norm(u)
    return([(u/r).tolist(),(u/-r).tolist()])

def kro(val1,val2):
    if val1==val2:
        return 1
    else:
        return 0

    
def bonval2(m,n):
    u = subtract(m,n).tolist()
    r = norm(u)
    u = (u/r).tolist()
    d2 = np.zeros((6,6))
    for i in range(2):
        for j in range(3):
            for k in range(2):
                for l in range(3):
                    tmp1=(u[j]*u[l]-kro(j,l))/r
                    if i==k:
                        tmp1*=-1.0
                    d2[3*i+j][3*k+l]=tmp1
    return d2
                    
# ---> Angle
def ifparll(list1,list2):
    plimit=1.0e-10
    if abs(abs(np.dot(list1,list2))-1.0e0) > plimit:
        return False
    else:
        return True

def dot_pro(list1,list2):
    return (list1[0]*list2[0]+list1[1]*list2[1]+list1[2]*list2[2])

def zeta(ww1,ww2,ww3):
    if ww1==ww2:
        return 1
    elif ww1==ww3:
        return -1
    else:
        return 0
        
def angval(m,o,n):
    u = subtract(m,o).tolist()
    v = subtract(n,o).tolist()
    r1 = norm(u)
    r2 = norm(v)
    u = (u/r1).tolist()
    v = (v/r2).tolist()
    # -- w vector --
    if not ifparll(u,v):
        wp=np.cross(u,v).tolist()
    else:
        if (not ifparll(u,[1,-1,1]) and not ifparll(v,[1,-1,1])):
            wp=np.cross(u,[1,-1,1])
        elif ( parll(u,[1,-1,1]) and parll(v,[1,-1,1])):
            wp=np.cross(u,[-1,1,1])
    
    wp = (wp/norm(wp)).tolist()
    uw = np.cross(u,wp)
    vw = np.cross(wp,v)
    dadx = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            dadx[i][j] = zeta(i,0,1)*uw[j]/r1 + zeta(i,2,1)*vw[j]/r2
    return dadx

def angval2(m,o,n):
    u = subtract(m,o).tolist()
    v = subtract(n,o).tolist()
    r1 = norm(u)
    r2 = norm(v)
    u = (u/r1).tolist()
    v = (v/r2).tolist()
    qa = acos(np.dot(u,v))
    # -- w vector --
    if not ifparll(u,v):
        wp=np.cross(u,v).tolist()
    else:
        if (not ifparll(u,[1,-1,1]) and not ifparll(v,[1,-1,1])):
            wp=np.cross(u,[1,-1,1])
        elif ( parll(u,[1,-1,1]) and parll(v,[1,-1,1])):
            wp=np.cross(u,[-1,1,1])
    
    wp = (wp/norm(wp)).tolist()
    uw = np.cross(u,wp).tolist()
    vw = np.cross(wp,v).tolist()
    #-- First derivative --
    dadx = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            dadx[i][j] = zeta(i,0,1)*uw[j]/r1 + zeta(i,2,1)*vw[j]/r2

    #-- Second derivative --
    da2dx2 = np.zeros((9,9))
    for a in range(3):
        for i in range(3):
            for b in range(3):
                for j in range(3):
                    tmp1 = (zeta(a,0,1)*zeta(b,0,1)*(u[i]*v[j] + u[j]*v[i] - \
                                                     3*u[i]*u[j]*cos(qa) + kro(i,j)*\
                                                     cos(qa)))/(r1**2)*sin(qa)
                    tmp1+= (zeta(a,2,1)*zeta(b,2,1)*(v[i]*u[j] + v[j]*u[i] - \
                                                     3*v[i]*v[j]*cos(qa) + kro(i,j)*\
                                                     cos(qa)))/(r2**2)*sin(qa)
                    tmp1+= (zeta(a,0,1)*zeta(b,2,1)*(u[i]*u[j] + v[j]*v[i] - u[i]*\
                                                     v[j]*cos(qa) - kro(i,j)))/r1*r2*sin(qa)
                    tmp1+= (zeta(a,2,1)*zeta(b,0,1)*(v[i]*v[j] + u[i]*u[j] - v[i]*\
                                                     u[j]*cos(qa) - kro(i,j)))/r1*r2*sin(qa)
                    tmp1-= (cos(qa)/sin(qa))*dadx[a][i]*dadx[b][j]
                    da2dx2[3*a+i][3*b+j] = tmp1
    return da2dx2

# ---> Torsional
    
def dihval(m,o,p,n):
    u = subtract(m,o).tolist()
    v = subtract(n,p).tolist()
    w = subtract(p,o).tolist()
    ru= norm(u)
    rv= norm(v)
    rw= norm(w)
    u = (u/ru).tolist()
    v = (v/rv).tolist()
    w = (w/rw).tolist()
    ucos = np.dot(u,w)
    vcos = -np.dot(v,w)
    usin = sqrt(1-(np.dot(u,w))**2)
    vsin = sqrt(1-(np.dot(v,w))**2)
    uw = np.cross(u,w).tolist()
    vw = np.cross(v,w).tolist()
    dddx = np.zeros((4,3))
    for a in range(4):
        for i in range(3):
            tmp1=tmp2=tmp3=tmp4=0.0
            if (a==0 or a==1):
                tmp1 = zeta(a,0,1)*uw[i] / (ru*usin*usin)
            if (a==2 or a ==3):
                tmp2 = zeta(a,2,3)*vw[i] / (rv*vsin*vsin)
            if (a==1 or a==2):
                tmp3 = zeta(a,1,2)*uw[i]*ucos / (rw * usin * usin)
            if (a==1 or a==2):
                tmp4 =-zeta(a,2,1)*vw[i]*vcos / (rw * vsin * vsin)
            dddx[a][i] = tmp1 + tmp2 + tmp3 + tmp4
    return dddx

def dihval2(m,o,p,n):
    u = subtract(m,o).tolist()
    v = subtract(n,p).tolist()
    w = subtract(p,o).tolist()
    ru= norm(u)
    rv= norm(v)
    rw= norm(w)
    u = (u/ru).tolist()
    v = (v/rv).tolist()
    w = (w/rw).tolist()
    ucos = np.dot(u,w)
    vcos = -np.dot(v,w)
    usin = sqrt(1-(np.dot(u,w))**2)
    vsin = sqrt(1-(np.dot(v,w))**2)
    usin4= usin*usin*usin*usin
    vsin4= vsin*vsin*vsin*vsin
    ucos3= ucos*ucos*ucos
    vcos3= vcos*vcos*vcos
    uw = np.cross(u,w).tolist()
    vw = np.cross(v,w).tolist()
    dd2dx2 = np.zeros((12,12))
    for a in range(4):
        for b in range(4):
            for i in range(3):
                for j in range(3):
                    tmp1=0
                    if ((a==0 and b==0) or (a==1 and b==0) or (a==1 and b==1)):
                        tmp1 += (zeta(a,0,1)*zeta(b,0,1)*(uw[j]*(w[j]*ucos-u[i]) + \
                                            uw[j]*(w[i]*ucos-u[i])))/(ru*ru*usin4)
                    if ((a==3 and b==3) or (a==3 and b==2) or (a==2 and b==2)):
                        tmp1 += (zeta(a,3,2)*zeta(b,3,2)*(vw[i]*(w[j]*vcos+v[j]) + vw[j]*\
                                                    (w[i]*vcos+v[i])))/(rv*rv*vsin4)
                    if ((a==1 and b==1) or (a==2 and b==1) or (a==2 and b==0) or (a==1 and b==0)):
                        tmp1 += ((zeta(a,0,1)*zeta(b,1,2)+zeta(a,2,1)*zeta(b,1,0))*(uw[i]*\
                                            (w[j]-2*u[j]*ucos+w[j]*ucos*ucos)+uw[j]*\
                                        (w[i]-2*u[i]*ucos+w[i]*ucos*ucos)))/(2*ru*rw*usin4)
                    if ((a==3 and b==2) or (a==3 and b==1) or (a==2 and b==2) or (a==2 and b==1)):
                        tmp1 += ((zeta(a,3,2)*zeta(b,2,1)+zeta(a,1,2)*zeta(b,2,3))*(vw[i]*\
                            (w[j]+2*v[j]*vcos+w[j]*vcos*vcos)+vw[j]*\
                            (w[i]+2*v[i]*vcos+w[i]*vcos*vcos)))/(2*rv*rw*vsin4)
                    if ((a==1 and b==1) or (a==2 and b==2) or (a==2 and b==1)):
                        tmp1 += (zeta(a,1,2)*zeta(b,2,1)*(uw[i]*(u[j]+u[j]*ucos*ucos-3*w[j]*ucos+w[j]*\
                            ucos3)+uw[j]*(u[i]+u[i]*ucos*ucos-3*w[i]*\
                            ucos+w[i]*ucos3)))/(a*rw*rw*usin4)
                    if ((a==2 and b==1) or (a==2 and b==2) or (a==1 and b==1)):
                        tmp1 += (zeta(a,2,1)*zeta(b,1,2)*(vw[i]*(-v[j]-v[j]*vcos*vcos-3*w[j]*\
                            vcos+w[j]*vcos3)+vw[j]*(-v[i]-v[i]*vcos*vcos-3*\
                        w[i]*vcos+w[i]*vcos3)))/(2*rw*rw*vsin4)
                    if (not a==b and not i==j):
                        if ( not i==0 and not j==0):
                            k=0
                        elif (not i==1 and not j==1):
                            k=1
                        else:
                            k=2
                        if (a==1 and b==1):
                            tmp1 += (zeta(a,0,1)*zeta(b,1,2)*(j-1)*pow(-0.5,fabs(j-i))*\
                                     (+w[k]*ucos-u[k]))/(ru*rw*usin*usin)
                        if ((a==3 and b==2) or (a==3 and b==1) or (a==2 and b==2) or (a==2 and b==1)):
                            tmp1 += (zeta(a,3,2)*zeta(b,2,1)*(j-1)*pow(-0.5,fabs(j-i))*\
                                     (-w[k]*vcos-v[k]))/(rv*rw*vsin*vsin)
                        if ((a==2 and b==1) or (a==2 and b==0) or (a==1 and b==1) or (a==1 and b==0)):
                            tmp1 += (zeta(a,2,1)*zeta(b,1,0)*(j-i)*pow(-0.5,fabs(j-i))*\
                                     (-w[k]*ucos+u[k]))/(ru*rw*usin*usin)
                        if (a==2 and b==2):
                            tmp1 += (zeta(a,1,2)*zeta(b,2,3)*(j-i)*pow(-0.5,fabs(j-i))*\
                                     (+w[k]*vcos+v[k]))/(rv*rw*vsin*vsin)
                    dd2dx2[3*a+i][3*b+j] = dd2dx2[3*b+j][3*a+i] = tmp1
    return dd2dx2
