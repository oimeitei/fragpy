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
Simple stuffs
'''

from numpy.linalg import norm
from numpy import subtract
from math import *
from fragpy.opt_param import param
from fragpy.odict import Ang
from copy import deepcopy
import numpy as np

def rmsB(list1):
    rmst = 0
    greatt = abs(list1[0])
    for i in list1: 
        if abs(i)>greatt:
            greatt = abs(i)
        rmst += i**2
    rmst = rmst/len(list1)
    rmst = rmst**0.5
    return rmst

def maxB(list1):
    greatt = abs(list1[0])
    for i in list1:
        if abs(i)> greatt:
            greatt = abs(i)
    return greatt

def measure_cart(list1):
    grad=0
    maxi = abs(list1[0][0])
    for i in list1:
        for j in range(3):
            grad += i[j]**2
            if abs(i[j]) > maxi:
                maxi = abs(i[j])
    grad /= param.natom*3
    grad = grad**0.5
    return(grad,maxi)

def measure_flat(L1):
    grad = 0.0
    maxi = abs(L1[0])
    for i in L1:
        grad += i**2
        if abs(i) > maxi:
            maxi = abs(i)
    grad /= param.natom*3
    grad = grad**0.5
    return(grad,maxi)

def sub_cart(list1,list2):
    val = len(list1)
    empty = np.zeros((val,3))
    for i in range(val):
        empty[i][0] = abs(abs(list1[i].qx) - abs(list2[i].qx))
        empty[i][1] = abs(abs(list1[i].qy) - abs(list2[i].qy))
        empty[i][2] = abs(abs(list1[i].qz) - abs(list2[i].qz))
    return empty

def makeblist(nf):
     bondlist = []
     bond_tmp = []
     bcoef = []
     
     for i in range(nf-1):
          for j in bondlist:
               if i+1 in j:
                    bond_tmp.append(i+1)
               elif i+2 in j:
                    bon_tmp.append(i+2)
          bondlist.append([i+1,i+2])
          bcoef.append(1)
     for i in bond_tmp:
          bondlist.append([i,0])
          bcoef.append(-1)
     return(bondlist,bcoef)

def makenblist(nf):
     # ---> Non-Bonded fragments
     nbondlist = []

     for i in range(nf):
          for j in range(i,nf):
               if not j+3 > nf:
                    nbondlist.append([i+1,j+3])
     return nbondlist


def neglect_nb(nbondlist,Thres,suplist):
     # ---> Neglect frag distant non-bonded fragments
     nbon_tmp = []
     for i1,i2 in nbondlist:
          min1 = 1000.0
          tmp1 = []; tmp2 = []
          for i in suplist:
               if i.grp == i1: tmp1.append(i)
               elif i.grp == i2: tmp2.append(i)
          for i in tmp1:
               for j in tmp2:
                    d = norm(subtract(i.q,j.q).tolist())
                    if d < min1:
                         min1 = d
          if min1 < Thres/Ang:
               nbon_tmp.append([i1,i2])
     nbondlist = deepcopy(nbon_tmp)
     return nbondlist
