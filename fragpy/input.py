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
Ang = 0.52918

def input_scan(inp):
    r=open(inp,'r')
    geom = False
    bohr = False
    cart_coord = []
    c1 = 1.0
    check1 = False
    check2 = []
    for line in r.readlines():

        if 'geometry' in line and not '!' in line:
            geom = True
            if line.split()[1] == 'Bohr':
                bohr = True
            continue

        if geom:
            if line.isspace():
                geom=False
                continue
            s=line.split()
            if bohr:
                c1 = Ang
            if len(s) == 4:
                cart_coord.append([s[0],c1*float(s[1]),c1*float(s[2]),c1*float(s[3])])
            elif len(s) == 5:
                cart_coord.append([s[0],c1*float(s[1]),c1*float(s[2]),c1*float(s[3]),int(s[4])])

        if 'set' in line and not '!' in line:
            check1 = True
            continue
        elif 'optf' in line and not '!' in line:
            check2.append('OPTF')

        if check1:
            if line.isspace():
                check1 = False
                continue
            
            s1=line.split('=')
            s1[1] = s1[1].replace('\n','')
            if s1[0] == 'con':          param.con = s1[1]                     
            if s1[0] == 'DEBUG':        param.DEBUG = int(s1[1])                           
            if s1[0] == 'Prestart':     param.Prestart =  True if s1[1]=='True' else False        
            if s1[0] == 'restart':      param.restart =  True if s1[1]=='True' else False        
            if s1[0] == 'basis':        param.basis =  s1[1]                   
            if s1[0] == 'basis_jk':     param.basis_jk = s1[1]           
            if s1[0] == 'basis_mp2':    param.basis_mp2 = s1[1]
            if s1[0] == 'cbasis':       param.cbasis = s1[1]
            if s1[0] == 'cbasis_jk':    param.cbasis_jk = s1[1]
            if s1[0] == 'cbasis_mp2':   param.cbasis_mp2 = s1[1]
            if s1[0] == 'E_conv':       param.E_conv = s1[1]
            if s1[0] == 'D_conv':       param.D_conv = s1[1]
            if s1[0] == 'O_conv':       param.O_conv = s1[1]
            if s1[0] == 'MAXG':         param.MAXG = float(s1[1])
            if s1[0] == 'DELE':         param.DELE = float(s1[1])
            if s1[0] == 'MAXQ':         param.MAXQ = float(s1[1])
            if s1[0] == 'bonthres':     param.bonthres = float(s1[1])
            if s1[0] == 'cbonthres':    param.cbonthres = float(s1[1])
            if s1[0] == 'thresnb':      param.thresnb  = float(s1[1])
            if s1[0] == 'cthresnb':     param.cthresnb = float(s1[1])
            if s1[0] == 'old_trust':    param.old_trust = True if s1[1]=='True' else False
            if s1[0] == 'OPTstep':      param.OPTstep= s1[1]
            if s1[0] == 'modelHtype':   param.modelHtype = s1[1] 
            if s1[0] == 'LBFGS':        param.LBFGS = True if s1[1]=='True' else False
            if s1[0] == 'ITMAX':        param.ITMAX=int(s1[1])     
            if s1[0] == 'trust':        param.trust = float(s1[1])
            if s1[0] == 'stepmax':      param.stepmax = float(s1[1])
            if s1[0] == 'stepmin':      param.stepmin = float(s1[1])
            if s1[0] == 'updateT':      param.updateT = True if s1[1]=='True' else False
            if s1[0] == 'updateFT':     param.updateFT =  True if s1[1]=='True' else False
            if s1[0] == 'updateH':      param.updateH = True if s1[1]=='True' else False
            if s1[0] == 'FIXGEOM':      param.FIXGEOM = True if s1[1]=='True' else False
            if s1[0] == 'l1o1':         param.l1o1 = True if s1[1]=='True' else False
            if s1[0] == 'nbscale':      param.nbscale = float(s1[1])
            if s1[0] == 'coarse':       param.coarse = True if s1[1]=='True' else False
            if s1[0] == 'coarse_level':  param.coarse_level =int(s1[1])     
            if s1[0] == 'dumfrag':      param.dumfrag = True if s1[1]=='True' else False
            if s1[0] == 'dumfrag_level':param.dumfrag_level =int(s1[1]) 
            if s1[0] == 'parallel':     param.parallel = int(s1[1]) 
            if s1[0] == 'procs':        param.procs = int(s1[1])   
            if s1[0] == 'program':      param.program = s1[1]
            if s1[0] == 'memory':       param.memory = int(s1[1])
            if s1[0] == 'MOLPRO':       param.MOLPRO = s1[1]
            if s1[0] == 'PSI4':         param.PSI4 = s1[1]
            if s1[0] == 'method':       param.method = s1[1]
            if s1[0] == 'cmethod':      param.cmethod = s1[1]
            if s1[0] == 'GAUSSIAN':     param.GAUSSIAN = s1[1]
            if s1[0] == 'gau_command':  
                s2 = s1[1].split(',')
                for comm in s2:
                    param.gau_command.append(comm)
    r.close()
    return(cart_coord,check2)
