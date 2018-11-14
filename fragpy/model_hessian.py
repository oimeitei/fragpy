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
import numpy as np 
from fragpy.opt_param import param
from fragpy.odict import *
from numpy.linalg import norm
from numpy import subtract

"""\
   Simple diagonal model Hessian
"""

def modelH(Q,Int):

    intcoord = Int.intcoord
    tmpclist = Int.tmpclist
    anglist  = Int.anglist 
    dih      = Int.dih     
    hbonlist = Int.hbonlist
    Hangle   = Int.Hangle  
    Hdih     = Int.Hdih     
    
 
    if param.modelHtype == 'R':
        mH = []
        # R. Lindh k(r) = 0.45
        for i in tmpclist:
            Q1 = Q[int(i[0])-1]; Q2 = Q[int(i[1])-1]
            mH.append(0.45*RHO(Q1.q,Q2.q,Q1.atm,Q2.atm))

        for i in anglist:
            Q1 = Q[int(i[0])-1]; Q2 = Q[int(i[1])-1]; Q3 = Q[int(i[2])-1] 
            mH.append(0.15 * RHO(Q1.q,Q2.q,Q1.atm,Q2.atm) * RHO(Q2.q,Q3.q,Q2.atm,Q3.atm))

        for i in dih:
            Q1 = Q[int(i[0])-1]; Q2 = Q[int(i[1])-1]; Q3 = Q[int(i[2])-1]; Q4 = Q[int(i[3])-1]
            mH.append(0.005 * RHO(Q1.q,Q2.q,Q1.atm,Q2.atm) * RHO(Q2.q,Q3.q,Q2.atm,Q3.atm) *
                      RHO(Q3.q,Q4.q,Q3.atm,Q4.atm))

        for i in hbonlist:
            Q1 = Q[int(i[0])-1]; Q2 = Q[int(i[1])-1]
            mH.append(0.45*RHO(Q1.q,Q2.q,Q1.atm,Q2.atm))

        for i in Hangle:
            Q1 = Q[int(i[0])-1]; Q2 = Q[int(i[1])-1]; Q3 = Q[int(i[2])-1] 
            mH.append(0.15 * RHO(Q1.q,Q2.q,Q1.atm,Q2.atm) * RHO(Q2.q,Q3.q,Q2.atm,Q3.atm))

        for i in Hdih:
            Q1 = Q[int(i[0])-1]; Q2 = Q[int(i[1])-1]; Q3 = Q[int(i[2])-1]; Q4 = Q[int(i[3])-1]
            mH.append(0.005 * RHO(Q1.q,Q2.q,Q1.atm,Q2.atm) * RHO(Q2.q,Q3.q,Q2.atm,Q3.atm) *
                      RHO(Q3.q,Q4.q,Q3.atm,Q4.atm))

        tempH = np.eye(intcoord,dtype=np.float64)
        tempH = mH*tempH
    elif param.modelHtype == 'S':
        # Theoret. Chim. Acta 66, 333, 1984
        mH = []
        for i in tmpclist:
            Q1 = Q[int(i[0])-1]; Q2 = Q[int(i[1])-1]
            mH.append(Fstr(Q1.q,Q2.q,Q1.atm,Q2.atm))    

        for i in anglist:
            Q1 = Q[int(i[0])-1]; Q3 = Q[int(i[2])-1] 
            mH.append(Fbend(Q1.atm,Q3.atm))

        for i in dih:
            Q2 = Q[int(i[1])-1]; Q3 = Q[int(i[2])-1]
            mH.append(Ftors(Q2.q,Q3.q,Q2.atm,Q3.atm))

        for i in hbonlist:
            Q1 = Q[int(i[0])-1]; Q2 = Q[int(i[1])-1]
            mH.append(Fstr(Q1.q,Q2.q,Q1.atm,Q2.atm))

        for i in Hangle:
            Q1 = Q[int(i[0])-1]; Q3 = Q[int(i[2])-1] 
            mH.append(Fbend(Q1.atm,Q3.atm))
        for i in Hdih:
            Q2 = Q[int(i[1])-1]; Q3 = Q[int(i[2])-1]
            mH.append(Ftors(Q2.q,Q3.q,Q2.atm,Q3.atm))

        tempH = np.eye(intcoord,dtype=np.float64)
        tempH = mH*tempH      
    else:
        out_error(' ERROR: Set modelHtype to eiter \'R\' or \'S\'\n')
    if param.DEBUG >= 1:
        out = open('output','a')
        out.write('\n')
        out.write(' Model Hessian initiated.\n')
        out.write('\n')
        out.close()
    return tempH.tolist()

    
def Fstr(la1,la2,sy1,sy2):
    F_r = norm(subtract(la1,la2))
    Fstr_1 = 1.734 / (F_r - getB(sy1,sy2))**3
    return Fstr_1

def Fbend(sy1,sy2):
    if sy1 == 'H' or sy2 == 'H':
        Fbend_1 = 0.160
    if (not sy1 == 'H' and not sy2 == 'H'):
        Fbend_1 = 0.250
    return Fbend_1

def Ftors(la1,la2,sy1,sy2):
    F_r = norm(subtract(la1,la2))
    C_r = covrad[sy1]+covrad[sy2]
    value1 = 0.0023 - (0.07*(F_r - C_r))
    return value1


def getAlpha(atm1,atm2):
    alpha_t_IJ = [[1.000,0.3949,0.3949],[0.3949,0.2800,0.2800],[0.3949,0.2800,0.2800]]
    PER1 = {'H':0,'C':1,'N':1,'O':1}
    value1 = alpha_t_IJ[PER1[atm1]][PER1[atm2]]
    return value1

def getr2(atm1,atm2):
    alpha_t_IJ = [[1.35,2.10,2.53],[2.10,2.87,3.40],[2.53,3.40,3.40]]
    PER1 = {'H':0,'C':1,'N':1,'O':1}
    value1 = alpha_t_IJ[PER1[atm1]][PER1[atm2]]
    return value1
    
def RHO(la1,la2,sy1,sy2):
    # IJQC, 106, 2536, 2006
    #np.exp(-((norm(subtract(la1,la2))/(covrad[sy1]+covrad[sy2]))-1))
    # R. Lindh Chem. Phys. Lett., 241, 423, 1995
    rhoIJ = np.exp(getAlpha(sy1,sy2)*(getr2(sy1,sy2)**2-(norm(subtract(la1,la2)))**2))
    return rhoIJ

def getB(atm1,atm2):
    B_array = [[-0.244,0.352],[0.352,1.085]]
    PER1 = {'H':0,'C':1,'N':1,'O':1}
    value1 = B_array[PER1[atm1]][PER1[atm2]]
    return value1
