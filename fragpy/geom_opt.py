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
import re
from fragpy.opt_param import param
from fragpy.bon_angle_dih import *
from fragpy.odict import *
from math import *

'''
Generate redundant internal coordinates
Bond, Hydrogen Bond, Angle & Torsional
'''
checker1=[]
checker2=[]
checker3=0

def internal_coord(Q):
    global nhlist
    nhlist, hlist = sorth(Q)
    tmpclist = connect(nhlist)
    hbonlist, lenHbon = hbon(hlist,nhlist)
    Hangle = hangle(hbonlist,tmpclist,lenHbon)
    anglist = angl(tmpclist)
    Hdih = hdih(Hangle,anglist)
    dih = dihed(anglist)
    intcoord = len(tmpclist)+len(anglist)+len(dih)+len(hbonlist)+len(Hangle)+len(Hdih)
    if param.DEBUG >= 1:
        out = open('output','a')
        out.write('-----------------------------\n')
        out.write('| No. of Internal Coordinates \n')
        out.write('| Bonds     : {:>5}\n'.format(len(tmpclist)+len(hbonlist)))
        out.write('| H bonds   : {:>5}\n'.format(lenHbon))
        out.write('| Angles    : {:>5}\n'.format(len(anglist)+len(Hangle)))
        out.write('| Dihedrals : {:>5}\n'.format(len(dih)+len(Hdih)))
        out.write('| Total     : {:>5}\n'.format(intcoord))
        out.write('-----------------------------\n')
        out.write('\n')
        out.close()
    return INTCOORD(intcoord,tmpclist,anglist,dih,hbonlist,Hangle,Hdih)

class INTCOORD:
    
    def __init__(self,intcoord,tmpclist,anglist,dih,hbonlist,Hangle,Hdih):
        self.intcoord = intcoord
        self.tmpclist = tmpclist
        self.anglist = anglist
        self.dih = dih
        self.hbonlist = hbonlist
        self.Hangle = Hangle
        self.Hdih = Hdih
        
def sorth(list1):
    list2=[]
    list3=[]
    for i in list1:
        if i.atm =='H':
            list2.append(i)
        else:
            list3.append(i)
    return(list3,list2)


def connect(l1): 
    empty = []
    for i in l1:
        for j in l1:
            if not i == j:
                d1 = GetBond(i.q,j.q)
                d2 = covrad[i.atm]+covrad[j.atm]
                if d1 <= d2+0.1:
                    if (not [i.n,j.n] in empty and not [j.n,i.n] in empty):
                        empty.append([i.n,j.n])
    return empty


def readfile(F1,coarse):
    suplist = []
    cout = 1
    nfrag = 3
    cnfrag = 1
    for i in F1:
        xq = i[1]/Ang; yq = i[2]/Ang; zq = i[3]/Ang
        grp = int(re.findall(r'\d+',i[0])[0]); atm = i[0][0]
        if coarse:
            cgrp = int(i[4])
        else:
            cgrp = 1
        if grp > nfrag:
            nfrag = grp
        if cgrp > cnfrag:
            cnfrag = cgrp            
        suplist.append(Coords(cout,[xq,yq,zq],atm,grp,cgrp))
        cout += 1
    return(suplist,nfrag,cnfrag)

class Coords:
    def __init__(self,n,q,atm,grp,cgrp):
        self.n = n
        self.qx = q[0]
        self.qy = q[1]
        self.qz = q[2]
        self.atm = atm
        self.q = [self.qx,self.qy,self.qz]
        self.grp = grp
        self.cgrp = cgrp
