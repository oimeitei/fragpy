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
Prints the output file
'''

import pickle
from fragpy.odict import Ang

def out_head(out,con):
    out.write('------------------------------------------------------------------\n')
    out.write('------------------------------------------------------------------\n')
    out.write('--    PROGRAM :: OPTF BY O. R. MEITEI                           --\n')
    out.write('--    Geometry optimizer using molecular fragments              --\n')
    out.write('------------------------------------------------------------------\n')
    out.write('\n')
    out.write('\n')
    if con=='s':
        out.write('                     Supermolecular optimization                  \n')
        out.write('\n')
    elif con=='fnb':
        out.write(' Fragment optimization with Bonding and Non-Bonding contribution  \n')
        out.write('             with Incremental correction at highest level         \n')
        out.write('Level of fragmentation : Level-1\n')
        out.write('Incremental correction : HF\n')
        out.write('\n')
    elif con=='fb':
        out.write('      Fragment optimization with only the Bonding contribution \n')
        out.write('             with Incremental correction at highest level         \n')
        out.write('Level of fragmentation : Level-1\n')
        out.write('Incremental correction : HF\n')
        out.write('\n')

def out_step(n):
    out=open('output','a')
    out.write('---------------------------\n')
    out.write('| Optimization step : {:3} |\n'.format(n))
    out.write('---------------------------\n')
    out.write('\n')
    out.write(' Starting Fragment calculations\n')
    out.close()

def out_opt(etot,diff,rmsg,maxg,rmsq,maxq,xrmsg,xmaxg,xrmsX,xmaxX):
    out=open('output','a')
    out.write('-------------------------------------------'\
              '-------------------------------------\n')
    out.write(' Energy          Del(E)     RMS(grad) '\
              ' MAX(grad)  RMS(step)  MAX(step)    \n')
    out.write('-------------------------------------------'\
              '-------------------------------------\n')
    out.write(' {:>13.8f}  {:>.3e}  {:>.3e}  {:>.3e}  '
              '{:>.3e}  {:>.3e} {:>4}\n'.format(etot,\
                              diff,rmsg,maxg,rmsq,maxq,'int'))
    out.write('                            {:>.3e}  {:>.3e}  '\
              '{:>.3e}  {:>.3e} {:>4}\n'.format(xrmsg,xmaxg,xrmsX,\
                                                xmaxX,'cart'))
    out.write('-------------------------------------------'\
              '-------------------------------------\n')
    out.close()

def out_stop(n):
    out = open('output','a')
    out.write(' Optmization succesful\n')
    out.write(' Optimization converged in {:3i} steps\n'.format(n))
    out.close()

def out_noconv():
    out = open('output','a')
    out.write('----------------\n')
    out.write(' NO CONVERGENCE\n')
    out.write('----------------\n')
    out.write('\n')
    out.close()

def out_error(*args):
    out = open('output','a')
    out.write('\n')
    for i in args:
        out.write(i)
    out.write('\n')
    out.close()
    quit()

class Store:
    def __init__(self):
        self.step = []
    
    def addS(self,s,e,de,gr,gm,sr,sm):
        self.step.append(SMethod(s,e,de,gr,gm,sr,sm))

    def printS(self):
        out = open('output','a')
        out.write(' Printing in a.u\n')
        out.write('\n')
        out.write('-----------------------------------------------------'\
                  '------------------------------\n')
        out.write(' Step   Energy         del(E)      RMS(grad)   MAX(grad) '\
                  '  RMS(step)   MAX(step)\n')
        out.write('------------------------------------------------------'\
                  '-----------------------------\n')
        for i in self.step:
            out.write(' {:>2}  {:>13.8f}  {:>12.8f}  {:>10.8f}  {:>10.8f} '\
                      ' {:>10.8f}  {:>10.8f}\n'.format(i.step,i.E,i.dE,i.grms,i.gmax,
                                           i.srms,i.smax))
        out.write('--------------------------------------------------------'\
                  '---------------------------\n')
        out.write('\n')
        out.close()

class SMethod:
    def __init__(self,s,e,de,gr,gm,sr,sm):
        self.step = s; self.E = e; self.dE = de; self.grms = gr
        self.gmax = gm; self.srms = sr; self.smax = sm

def restartin(f1):
    with open('restart.file','wb') as f:
        pickle.dump(f1,f)

def restartout():
    with open('restart.file','rb') as f:
        f1 = pickle.load(f)
    return f1

class Rstore:
    def __init__(self,stp,slist,int1,Qs,Gs,trst,hess,eprev,store):
        self.Ostep = stp; self.suplist = slist; self.int_coords = int1
        self.Qstore = Qs; self.Gstore = Gs; self.trust = trst
        self.hessianQ = hess; self.eprev = eprev; self.store = store

def out_fgeom(Q):
    out = open('output','a')
    out.write(' Final Geometry in Angstrom\n')
    out.write('\n')
    for i in Q:
        out.write('{:>3}   {:>10.6f}   {:>10.6f}   {:>10.6f}\n'.
                  format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
    out.write('\n')
    out.write(' And its\'s done !!! \n')
    out.write('\n') 
    out.close()

def out_log(Q,step):
    out = open('logfile','a')
    out.write(str(len(Q))+'\n')
    out.write('Step : '+str(step)+'\n')
    for i in Q:
        out.write('{:>3}   {:>10.6f}   {:>10.6f}   {:>10.6f}\n'.
                  format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
    out.write('\n')
    out.close()
