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
Geometry optimiser with RFO
Simple & Automated step restrictions
'''

import numpy.linalg as linalg
from numpy.linalg import norm
from fragpy.opt_param import param
from numpy import *
from math import *
import numpy as np
import numpy

def minimize(Qglist,hessianQ,n,trust):
    if param.OPTstep == 'RF':
        xiQ, delE = RF(Qglist,hessianQ,n,trust)
    elif param.OPTstep == 'SRF':
        xiQ, delE = SRF(Qglist,hessianQ,n,trust)
    else:
        out_error('ERROR: set OPTstep to either RF or SRF\n')
    return(xiQ,delE)


def RF(Qglist,hessianQ,n,trust):
    
    # Augmented hessian
    AH = np.zeros((n+1,n+1))
    for i in range(n):
        for j in range(n):
            AH[i][j] = hessianQ[i][j]
        AH[i][n] = AH[n][i] = Qglist[i]  
    
    # Lowest eigenvec, scale the last element to 1
    egval, evec = linalg.eigh(AH)
    QMIN = evec[:,0]
    for i in range(n):
        QMIN[i] /= QMIN[n]
    xiQ = [i for i in QMIN[:-1]]
    
    # Simple step restriction using norm
    NORM = norm(xiQ)
    
    if NORM > trust:
        for i in range(n):
            xiQ[i] *= (trust/NORM)
    
    # Project E change
    g = np.dot(Qglist,xiQ)
    h = 0.5*(np.dot(np.dot(xiQ,hessianQ),xiQ))
    d = 1 + np.dot(xiQ,xiQ)
    delE = (g+h)/d

    if param.DEBUG >= 1:
        out = open('output','a')
        out.write('\n')
        out.write(' Taking RFO step with simple step restriction\n')
        out.write(' Norm of the unscaled RFO step : {:>.3e}\n'.format(NORM))
        out.write(' Trust radius for the opt step : {:>.3e}\n'.format(trust))
        out.write(' Projected energy change       : {:>.3e}\n'.format(delE))
        out.write('\n')
        out.close()              

    return(xiQ,delE)

def SRF(Qglist,hessianQ,n,trust):

    # Eigenpair of H
    evalH, evecH = numpy.linalg.eigh(hessianQ)

    # Augmented hessian
    AH = np.zeros((n+1,n+1))
    for i in range(n):
        for j in range(n):
            AH[i][j] = hessianQ[i][j]
        AH[i][n] = AH[n][i] = Qglist[i]  

    alpha = 1
    SAH = np.zeros((n+1,n+1))
    Aco  = False
    if param.DEBUG >= 1:
        out = open('output','a') 
        out.write(' Taking Automated step restriction RFO step \n')
        out.write('\n')
        out.write(' Trust radius : '+str(trust)+'\n') 
    for it in range(21):
        if it == 20:
            if not sqrt(xiQxiQ) < param.stepmax:
                Aco = True
            break
        for i in range(n+1):
            for j in range(n):
                SAH[j][i] = AH[j][i] / alpha
            SAH[n][i] = AH[n][i]
            
        # Eigen part of scaled aug-Hessian
        EVAL, EVEC = numpy.linalg.eig(SAH)
        idx = EVAL.argsort()
        VAL_MIN = EVAL[idx]
        QMIN = EVEC[:,idx]
        QMIN = QMIN[:,0].real
        xiQ = [ 0 for i in range(n) ]

        # Scale the last component to One
        for i in range(n):
            xiQ[i] += QMIN[i]/QMIN[n]
        xiQ = np.array(xiQ)
        
        if it == 0:
            savexiQ = xiQ
        
        # dq^tdq
        xiQxiQ = np.dot(xiQ,xiQ)

        # Analytical derivative d(dq^tdq)/d(alpha)
        lamda = np.dot(Qglist,xiQ)
        SUM = 0
        for i in range(n):
            SUM += (pow(np.dot(evecH[:,i].T,Qglist),2))/ \
                  (pow((evalH[i]-lamda*alpha),3))
        annalytic = 2*((lamda)/(1+alpha*xiQxiQ))*SUM
        alpha += 2* (((trust*sqrt(xiQxiQ))-xiQxiQ)/annalytic)
        if param.DEBUG >= 1:
            out.write(' {:>2}   {:5.3f}  {:>1.2e}\n'.format(it,sqrt(xiQxiQ),alpha)) #t

        if alpha > 1e8:
            # Case of alpha explode
            if (not sqrt(xiQxiQ) < param.stepmax and \
                not sqrt(np.dot(savexiQ,savexiQ)) < param.stepmax) :
                Aco = True
            if sqrt(xiQxiQ) < param.stepmax:
                xiQ = xiQ
            elif sqrt(np.dot(savexiq,savexiQ)) < param.stepmax:
                xiQ = savexiQ
            break
        elif sqrt(xiQxiQ) < (trust + 1e-5):
            break
    if Aco:
        xiQ = RF(Qglist,hessianQ,n,trust)

    # Project E change
    g = np.dot(Qglist,xiQ)
    h = 0.5*(np.dot(np.dot(xiQ,hessianQ),xiQ))
    d = 1 + np.dot(xiQ,xiQ)
    delE = (g+h)/d

    if param.DEBUG >= 1:
        out = open('output','a')
        out.write('\n')
        out.write(' Projected energy change       : {:>.3e}\n'.format(delE))
        out.write('\n')
        out.close()   

    return(xiQ,delE)
