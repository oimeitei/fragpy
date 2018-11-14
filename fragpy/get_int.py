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
Backtransformation to Cartesian coordinates
'''


from numpy import *
from copy import deepcopy
import numpy as np
from numpy import *
from fragpy.bon_angle_dih import *
from fragpy.simple_maths import rmsB
from fragpy.geom_opt import Coords
from fragpy.opt_param import param

def GetInt(glist,bmatRF,Gmat_inv,Pmat):
    grad = dot(Gmat_inv,dot(array(bmatRF),glist).T)
    gradq = dot(Pmat,grad)

    return gradq

def GetCart(list1,list2,list3,list4,int_coords):
    
    # back transformation of step (int -> cart)
    # C. Peng et. al. J. Comp. Chem.
    # xnew = xint + Binv * int_step

    intcoord = int_coords.intcoord 
    tmpclist = int_coords.tmpclist 
    anglist  = int_coords.anglist  
    dih      = int_coords.dih      
    hbonlist = int_coords.hbonlist 
    Hangle   = int_coords.Hangle   
    Hdih     = int_coords.Hdih     

    xflat = []
    for i in list2: xflat.extend(i.q)
    ilist = getintcoord(list2,int_coords) 
    itarget = []
    for it in range(intcoord):
         itarget.append(ilist[it]+list1[it])
    dqGC = deepcopy(list1)
    conli = False
    for it in range(20):
         dxnewGC = np.dot(list4,dqGC)
         xGC = np.zeros(param.natom*3)
         dxnewGC.tolist()
         for iz in range(param.natom*3):
              xGC[iz] = xflat[iz] + dxnewGC[iz]
         iGC = getQcoord(xGC,int_coords)  
         if it == 0:
              SAVEDX = deepcopy(xGC)
              SAVEDQ = deepcopy(iGC)
         rmsGC = rmsB(dxnewGC)
         
         if it == 0:
              prevGC = rmsGC
         diffGC = prevGC - rmsGC
         prevGC = rmsGC
         if diffGC < 0.0:
              xGC = SAVEDX
              iGC = SAVEDQ
              conli = True
              break
         if rmsGC < 1.0e-6:
              break
         dqGC = np.zeros(intcoord)
         # Dihedrals are transformed for comparison !
         cout = 0
         for iz in range(len(tmpclist)+len(anglist)):
              dqGC[cout] = (itarget[cout] - iGC[cout])
              cout += 1
         for iz in range(len(dih)):
              # transformation here
              tGC = dihtrans(itarget[cout],iGC[cout])   
              dqGC[cout] = (itarget[cout] - tGC)
              cout += 1
         for iz in range(len(hbonlist)+len(Hangle)):
              dqGC[cout] = (itarget[cout] - iGC[cout])
              cout += 1
         for iz in range(len(Hdih)):
              # transformation here
              tGC = dihtrans(itarget[cout],iGC[cout])    
              dqGC[cout] = (itarget[cout] - tGC)
              cout += 1 
         xflat = deepcopy(xGC)

    emptylist = []
    for it in range(param.natom): # convert to Coords object
         emptylist.append(Coords(it+1,[xGC[3*it+0],xGC[3*it+1],
                                       xGC[3*it+2]],list2[it].atm,
                                 list2[it].grp,list2[it].cgrp))
    return emptylist


def getintcoord(list1gic,int_coords):

    tmpclist = int_coords.tmpclist
    anglist  = int_coords.anglist 
    dih      = int_coords.dih     
    hbonlist = int_coords.hbonlist
    Hangle   = int_coords.Hangle  
    Hdih     = int_coords.Hdih    
    
    intQ = []
    for i in tmpclist:
        intQ.append(GetBond(list1gic[int(i[0])-1].q,list1gic[int(i[1])-1].q))
    for i in anglist:
        intQ.append(GetAngle(list1gic[int(i[0])-1].q,list1gic[int(i[1])-1].q,list1gic[int(i[2])-1].q))
    for i in dih:
        intQ.append(GetDihedral(list1gic[int(i[0])-1].q,list1gic[int(i[1])-1].q,list1gic[int(i[2])-1].q,\
                                list1gic[int(i[3]-1)].q))
    for i in hbonlist:
        intQ.append(GetBond(list1gic[int(i[0])-1].q,list1gic[int(i[1])-1].q))
    for i in Hangle:
        intQ.append(GetAngle(list1gic[int(i[0])-1].q,list1gic[int(i[1])-1].q,list1gic[int(i[2])-1].q))
    for i in Hdih:
        intQ.append(GetDihedral(list1gic[int(i[0])-1].q,list1gic[int(i[1])-1].q,list1gic[int(i[2])-1].q,\
                                list1gic[int(i[3]-1)].q))
    return intQ

    
def getQcoord(list1gqc,int_coords):
    # flat to standard
    emptylistgqc = []
    for i in range(param.natom):
        emptylistgqc.append(QCOORD([list1gqc[3*i+0],list1gqc[3*i+1],list1gqc[3*i+2]]))
    intQcoord = getintcoord(emptylistgqc,int_coords)
    return intQcoord

def dihtrans(val1,val2):
     # https://doi.org/10.1186/s12859-017-1834-2
     if abs(val2-val1) <= radians(180):
          val3 = val2
     elif val2 - val1 <= -radians(180):
          val3 = val2 + radians(360)
     elif val2 - val1 >= radians(180):
          val3 = val2 - radians(360)
     return val3

class QCOORD:
    def __init__(self,q):
        self.q = q
