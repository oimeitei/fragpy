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
"""\
 Bond lengths;  angles; torsionals
"""

from fragpy.odict import *
from math import *
from numpy import subtract
from numpy.linalg import norm
import numpy as np


checker1=[]
checker2=[]
checker3=0

def hbon(list1,list2):

    tmpii=[]
    electNEG = ['O','N','F','Cl','S']
    for i in list1:
        closeN = closest(i,list2,checker1,checker2,checker3)
        for j in list2:
            if not closeN == j:
                if j.atm in electNEG:
                    bond1 = GetBond(i.q,j.q)
                    bondc = covrad[j.atm] + covrad[i.atm]
                    bondv = vanrad[j.atm] + vanrad[i.atm]
                    if (bond1 > bondc and bond1 < 0.9000*(bondv)):
                        angle1 = degrees(GetAngle(closeN.q,i.q,j.q))
                        if angle1 > 90:
                            tmpii.append([i.n,j.n])

    lentmpii = len(tmpii)
            
    for i in list2:
        hclose1=closest(i,list1,checker1,checker2,checker3)
        if not hclose1==[]:
            tmpii.append([hclose1.n,i.n])
            hclose2=closest(i,list1,hclose1,checker2,checker3)
            if not hclose2==[]:
                tmpii.append([hclose2.n,i.n])
                hclose3=closest(i,list1,hclose1,hclose2,hclose2)
                if not hclose3==[]:
                    tmpii.append([hclose3.n,i.n])

    return(tmpii,lentmpii)
               
def hangle(list1,list2,lenHbon):  # Hbondlist, nonHbondlist
    hangl=[]
    for i in list1[lenHbon:]:
        tmpii=[]
        for j in list2:
            if i[1]==j[0]:
                tmpii.append(j[1])
            elif i[1]==j[1]:
                tmpii.append(j[0])
                #tmpii=j[1]
        for j in tmpii:
            if not [j,i[1],i[0]] in hangl:
                hangl.append([i[0],i[1],j])
        tmpii = []
        for j in list1[lenHbon:]:
            if (i[1] == j[1] and not i[0] == j[0]):
                tmpii.append(j[0])
        for j in tmpii:
            if not [j,i[1],i[0]] in hangl:
                hangl.append([i[0],i[1],j])
    return hangl       

def hdih(list1,list2):
    hdi=[]
    for i in list1:
        tmpii = []
        for j in list2:
            if (i[1] == j[0] and i[2] == j[1]):
                tmpii.append([i[0],i[1],i[2],j[2]])
            elif (i[1] == j[2] and i[2] == j[1]):
                tmpii.append([i[0],i[1],i[2],j[0]])
        for j in tmpii:
            if not [j[3],j[2],j[1],j[0]] in hdi:
                hdi.append(j)
        tmpii = []
        for j in list1:
            if (j[2] == i[1] and j[1] == i[2]):
                tmpii.append([i[0],i[1],i[2],j[0]])
        for j in tmpii:
            if not [j[3],j[2],j[1],j[0]] in hdi:
                hdi.append(j) 
    return hdi

def angl(list1):
    condlist=[]
    tmpang=[]
    ange=[]
    coun=0
    for i in list1:
        tmpii = []
        for j in list1:
            if (i[0] == j[0] and not i[1] == j[1]):
                tmpii.append([i[1],i[0],j[1]])
            elif i[0] == j[1]:
                tmpii.append([i[1],i[0],j[0]])
            elif i[1] == j[0]:
                tmpii.append([i[0],i[1],j[1]])
            elif (i[1] == j[1] and not i[0] == j[0]):
                tmpii.append([i[0],i[1],j[0]])
        for j in tmpii:
            if not sorted(j) in tmpang:
                ange.append(j)
                tmpang.append(sorted(j))
    return ange
        
def dihed(list1):
    condlist=[]
    dihlist=[]
    for i in list1:
        tmpii = []
        for j in list1:
            if (i[1]==j[0] and i[2]==j[1]):
                tmpii.append([i[0],i[1],i[2],j[2]])
            elif (i[1]==j[2] and i[2]==j[1]):
                tmpii.append([i[0],i[1],i[2],j[0]])
            elif (i[1]==j[2] and i[0]==j[1]):
                tmpii.append([i[2],i[1],i[0],j[0]])
            elif (i[1]==j[0] and i[0]==j[1]):
                tmpii.append([i[2],i[1],i[0],j[2]])
        for j in tmpii:
            if not sorted(j) in condlist:
                dihlist.append(j)
                condlist.append(sorted(j))
    return dihlist

def closest(x,list,check1,check3,check2):
    condition=0
    for i in list:
        if not i==check1:
            if not i==check2:
                if not i == check3:
                    d1 = GetBond(x.q,i.q)
                    d2 = covrad[x.atm]+covrad[i.atm]
                    if d1 <= d2+0.1:
                        close=i
                        condition=1
                        break
                    else:
                        close=[]
                        condition=2
    if condition==0:
        close=[]
    return (close)

def GetBond(list1,list2):
    bond1 = subtract(list1,list2).tolist()
    bond2 = norm(bond1)
    return bond2
def GetAngle(s1,s2,s3):
    import numpy

    u = subtract(s1,s2)
    v = subtract(s3,s2)
    ru = norm(u)
    rv = norm(v)
    u = (u/ru).tolist()
    v = (v/rv).tolist()
    angle1 = numpy.arccos(np.dot(u,v))
    return angle1


def GetDihedral(p0,p1,p2,p3):
    b0 = -1.0 * subtract(p1,p0)
    b1 = subtract(p2,p1)
    b2 = subtract(p3,p2)
    b1 /= norm(b1)
    v = subtract(b0,np.dot(b0,b1)*b1)
    w = subtract(b2,np.dot(b2,b1)*b1)
    x = np.dot(v,w)
    y = np.dot(np.cross(b1,v),w)
    return np.arctan2(y,x)
    
