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
#!/usr/bin/env python
import sys
from fragpy.opt_param import param
from fragpy.out_print import * 
from fragpy.geom_opt import internal_coord,readfile
from fragpy.model_hessian import modelH
from fragpy.bmat import BmatB,PmatP
from fragpy.update import Hupdate,Tupdate,FTupdate
from fragpy.opt_step import minimize
from fragpy.main import mainfunc
from fragpy.get_int import *
from fragpy.simple_maths import makeblist,makenblist,neglect_nb,\
    rmsB,maxB,measure_cart,measure_flat,sub_cart
from fragpy.opt_param import param
import numpy as np 
from numpy import *

'''
Runs the fragment geometry optimization.

'''
def optf(supfile):
    # ---> Initialisation
    suplist,nfrag,cnfrag = readfile(supfile,param.coarse)
    param.natom = len(suplist)
    param.suplist = suplist
    
    param.nfrag = nfrag    # No. of groups
    param.cnfrag = cnfrag
              
    # ---> Fragment list
    bondedlist,bcoef = makeblist(nfrag)
    param.bondedlist = bondedlist
    param.bcoef = bcoef
    if param.coarse:
        cbondedlist,cbcoef = makeblist(cnfrag)
        param.cbondedlist = cbondedlist
        param.cbcoef = cbcoef
        param.cnbondedlist = makenblist(cnfrag)
    param.nbondedlist = makenblist(nfrag)
    param.nbondedlist = neglect_nb(param.nbondedlist,param.thresnb,param.suplist)
    if param.coarse:
        param.cnbondedlist = neglect_nb(param.cnbondedlist,param.cthresnb,param.suplist)
    
    out = open('output','w')
    out_head(out,param.con)
    out.close()

    out = open('logfile','w')
    out.close()
    
    convergence = False
    
    if not param.restart:
        Qstore = []; Gstore = []; delE = 0.0; trust = param.trust
    
        # Generate redundant internal coords
        int_coords = internal_coord(suplist)  
    
        # Initial Hessian
        hessianQ = modelH(suplist,int_coords)
        
        store = Store(); eprev = 0.0; Ostep = 0
    
    for i in range(param.ITMAX):
    
        Ostep += 1
        if param.Prestart:
            if not param.restart:
                rest = Rstore(Ostep,suplist,int_coords,Qstore,Gstore,trust,
                              hessianQ,eprev,store)
                restartin(rest)
            elif param.restart:
                if param.DEBUG >= 1:
                    out = open('output','a')
                    out.write(' Restarting the fragment optimisation\n')
                    out.write('\n')
                    out.close()
                reobj = restartout()
                Ostep = reobj.Ostep; suplist = reobj.suplist; int_coords = reobj.int_coords  
                Qstore = reobj.Qstore; Gstore = reobj.Gstore; trust      = reobj.trust       
                hessianQ = reobj.hessianQ; eprev = reobj.eprev; store = reobj.store
                restart = False
    
        if i > 1: param.updateT = True
        if Ostep > 1: param.updateH = True
    
        # Gradients & Energy 
        out_log(suplist,Ostep)
        out_step(Ostep)
        glist,etot = mainfunc(suplist)
    
        # Generate B matrix, pseudo inverse & projection matrix
        bmatRF, binvRF = BmatB(suplist,int_coords)
        Pmat, Gmat_inv = PmatP(bmatRF,int_coords.intcoord)
    
        #store
        if Ostep == 1:
            eprev = etot
        delE = eprev - etot; eprev = etot
        Qstore.append(suplist); Gstore.append(glist)
    
        # Update trust
        if param.updateT:
            if param.updateFT:
                trust = FTupdate(imaxg,delE)
            else:
                trust = Tupdate(delE,delE_proj,trust,xiQ)
    
        # Cart ->> internal coord (gradients)
        glistQ = GetInt(glist,bmatRF,Gmat_inv,Pmat)
    
        # Update Hessian
        if param.updateH:
            hessianQ = Hupdate(hessianQ,Qstore,Gstore,int_coords,Pmat,suplist)
        else:
            # Project the hessian
            hessianQ = np.dot(Pmat,np.dot(hessianQ,Pmat))
        
        # Get Opt step
        xiQ, delE_proj = minimize(glistQ,hessianQ,int_coords.intcoord,trust)
    
        # Backtransform Opt step to Cart
        suplist = GetCart(xiQ,suplist,bmatRF,binvRF,int_coords)
    
        irmsg = rmsB(glistQ); imaxg = maxB(glistQ)
        irmsq = rmsB(xiQ); imaxq = maxB(xiQ)
        xiX = sub_cart(suplist,Qstore[-1])
        xrmsg, xmaxg = measure_flat(glist); xrmsq, xmaxq = measure_cart(xiX)
    
        out_opt(etot,delE,irmsg,imaxg,irmsq,imaxq,xrmsg,xmaxg,xrmsq,xmaxq)
        store.addS(Ostep,etot,delE,irmsg,imaxg,irmsq,imaxq)
    
        # New internal coordinates
        out = open('output','a'); out.write('\n')
        if not param.FIXGEOM:
            out.write(' New internal coordinates generated \n')
            int_coords = internal_coord(suplist)
        else:
            out.write(' No new internal coordinates generated \n')
        out.write('\n'); out.close()
    
        # Convergence check
    
        if ((imaxg < param.MAXG and param.DELE < 1.0e-6) or \
            (imaxg < param.MAXG and param.MAXQ < 1.8e-3)):
            convergence = True
            break
    
    if convergence:
        out_stop(Ostep)
    else:
        out_noconv()
    store.printS()
    out_fgeom(suplist)
