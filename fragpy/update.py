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
'''
Hessian and trust region update
'''

from fragpy.opt_param import param
from fragpy.bmat import BmatB,symm_matrix_inv
from fragpy.model_hessian import modelH
from fragpy.get_int import *
from math import *
from numpy.linalg import norm
import numpy as np


def Hupdate(hessianQ,Qstore,Gstore,int_coords,P_mat,suplist):
    
    if param.LBFGS:
        HessianQ = lbfgs(hessianQ,Qstore,Gstore,int_coords,P_mat,suplist)
        return HessianQ
    else:
        out_error('Other options not yet implemented for Hessian update, set LBFGS = True\n')
        

def lbfgs(hessianQ,Qstore,Gstore,int_coords,P_mat,suplist):

    intcoord = int_coords.intcoord
    tmpclist = int_coords.tmpclist
    anglist  = int_coords.anglist 
    dih      = int_coords.dih     
    hbonlist = int_coords.hbonlist
    Hangle   = int_coords.Hangle  
    Hdih     = int_coords.Hdih    
    
    # Empirical hessian
    hessianQ = modelH(suplist,int_coords)
    hessianQ = np.dot(P_mat,np.dot(hessianQ,P_mat))

    skipbfgs = 0

    if param.DEBUG >= 1:
        out = open('output','a')
        out.write('\n')
        out.write(' Limited memory BFGS update of the Hessian\n')
        out.write(' No. of steps used : {:4} Steps\n'.format(len(Gstore)))
        out.write('\n')
        out.close()
    for m in range(len(Gstore)-1): 

         # B matrix for 2 prev steps using current intcoord
         bmatRF1, binvRF1  = BmatB(Qstore[m+1],int_coords)
         bmatRF0, binvRF0  = BmatB(Qstore[m],int_coords)

         # G matrix & inverse
         G_mat1 = np.dot(np.array(bmatRF1),np.array(bmatRF1).T)
         G_inv1 = symm_matrix_inv(G_mat1,intcoord)
         G_mat0 = np.dot(np.array(bmatRF0),np.array(bmatRF0).T)
         G_inv0 = symm_matrix_inv(G_mat0,intcoord)

         # Gradients (projected) in the current int coords
         Gstore1 = GetInt(Gstore[m+1],bmatRF1,G_inv1,P_mat)
         Gstore0 = GetInt(Gstore[m],bmatRF0,G_inv0,P_mat)

         # Int coords in the current definition 
         Qstore1 = getintcoord(Qstore[m+1],int_coords)
         Qstore0 = getintcoord(Qstore[m],int_coords)
    
         xi = np.zeros(intcoord)
         dg = np.zeros(intcoord)

         for i in range(intcoord):
              dg[i] = Gstore1[i] - Gstore0[i]

         # dihedrals are transformed for comparison
         cout = 0
         for iz in range(len(tmpclist)+len(anglist)):
              xi[cout] = (Qstore1[cout] - Qstore0[cout])
              cout += 1
         for iz in range(len(dih)):
              # transformation here
              tGC = dihtrans(Qstore1[cout],Qstore0[cout])
              xi[cout] = (Qstore1[cout] - tGC)
              cout += 1
         for iz in range(len(hbonlist)+len(Hangle)):
              xi[cout] = (Qstore1[cout] - Qstore0[cout])
              cout += 1
         for iz in range(len(Hdih)):
              # transformation here
              tGC = dihtrans(Qstore1[cout],Qstore0[cout])
              xi[cout] = (Qstore1[cout] - tGC)
              cout += 1

         xi = np.array(xi)
         dg = np.array(dg)

         # Skip if Denominators dg*dq or dq*dq are less than 1.0e-7
         dgdq = np.dot(dg,xi)
         dqdq = np.dot(xi,xi)
         if dgdq < 1.0e-7 or dqdq < 1.0e-7 :
              skipbfgs += 1
              continue
         # Skip if max internal coord change > 0.5
         if max(xi) > 0.5:
              skipbfgs += 1
              continue
         DH = dg[None, :]*dg[:, None]/dot(xi,dg) - \
              hessianQ.dot(xi[None, :]*xi[:, None]).\
              dot(hessianQ)/xi.dot(hessianQ).dot(xi)
         # Restricting Hessian update to a limit
         for i in range(intcoord):
              for j in range(intcoord):
                  if fabs(0.5e0*DH[i][j]) > 1.0e0:
                       DH[i][j] = copysign(1.0,DH[i][j])   
         DH = np.dot(P_mat,np.dot(DH,P_mat)) #tmp
         hessianQ[:, :] = hessianQ + DH
    return hessianQ

# Trust update
def Tupdate(actual,projected,trust,xiQ):
    if not param.old_trust:
        Eratio = -actual/projected
        if Eratio < 0.25:
            if not trust/4.0 < param.stepmin:
                trust/= 4.00
            else:
                trust = param.stepmin
        elif Eratio > 0.75:
            if not trust*3.00 > param.stepmax:
                trust *= 3.00
            else:
                trust = param.stepmax
                
    else:
        Eratio = -actual/projected
        re1=0+0.75
        reu=2-0.75
        ri1=0+0.80
        riu=2-0.80
        if (Eratio < re1 or Eratio > reu):
            trust = (norm(xiQ)/2.0)
        elif (ri1 <=Eratio<=riu and abs(norm(xiQ)-trust)<=0.000188973):
            trust = trust*sqrt(2)
    if param.DEBUG:
        out = open('output','a')
        out.write('\n')
        out.write(' Trust updated to {:.3e} \n'.format(trust))
        out.write('\n')
        out.close()
    return trust

def FTupdate(imaxg,delE):

    # Fixed trust update based on 
    # maximum component of the gradient.
    # Exercise care with this 
    
    if param.DEBUG:
        out = open('output','a')
        out.write('\n')
        out.write(' Fixed trust update set ON !!! \n')
        out.write(' This is a hard fix on the trust depending on the max gradient\n')
        out.write('\n')

    if imaxg > 3.0e-3:
        dxmax = 0.3
    elif imaxg > 1.0e-3:
        dxmax = 0.2
    elif imaxg > 6.0e-4:
        dxmax = 0.1
    elif imaxg > 5.0e-4:
        dxmax = 0.05
    elif imaxg > 4.0e-4:
        dxmax = 0.025
    elif imaxg > 3.0e-4:
        dxmax = 0.01
    else:
        dxmax = 0.005

    # Further fix the trust if the energy rise
    if delE < 0.0e0:
        dxmax = 0.03
    
    return dxmax
