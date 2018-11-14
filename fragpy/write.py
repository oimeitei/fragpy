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
Creates input files for computing energies & gradients
using other codes
PSI4
Molpro

'''

from fragpy.opt_param import param
from fragpy.odict import *
from fragpy.out_print import out_error
from fragpy.geom_opt import Coords
from numpy import subtract
from numpy.linalg import norm
import numpy as np

'''
Creates fragments and writes input files
Also writes the geometry file 123.xyz
'''

def splitwritef(Q):

    if param.con == 's':
        if param.DEBUG >= 1:
            out = open('output','a')
            out.write(' Starting supermolecular optimisation '\
                      '- The conventional optimisation\n')
        if param.program.upper() == 'PSI4':
            PSIinp('psi-s.in',Q,'supMP2')
            return(['psi-s.in'],[])
        elif param.program.upper() == 'MOLPRO':
            MOLinp('mol-s.in',Q,'supMP2')
            return(['mol-s.in'],[])
        elif param.program.upper() == 'GAUSSIAN':
            GAUinp('gau-s.in',Q,'supMP2')
            return(['gau-s.in'],[])
        else:
            out_error('ERROR: No program to run Fragpy - set param.program in opt_param.py\n') 
    frags = Fragments()
    frags.addcoord(Q)
    frags.addcaps(Q)
    frags.adddums(Q)   
    if param.con == 'fnb':
        frags.addcoord_nb(Q)
        frags.addcaps_nb(Q)

    # ---> Sup for correction term
    if not param.coarse:
        if param.program.upper() == 'PSI4':
            PSIinp('psi-s-hf.in',Q,'supHF')
        elif param.program.upper() == 'MOLPRO':
            MOLinp('mol-s-hf.in',Q,'supHF')
        elif param.program.upper() == 'GAUSSIAN':
            GAUinp('gau-s-hf.com',Q,'supHF')  
        else:
            out_error('ERROR: No program to run Fragpy - set param.program in opt_param.py\n')  
    else:
        if param.coarse_level == 2:
            if param.program.upper() == 'PSI4':
                PSIinp('cpsi-s-hf.in',Q,'supHFcca')
            elif param.program.upper() == 'MOLPRO':
                MOLinp('cmol-s-hf.in',Q,'supHFcca')
            elif param.program.upper() == 'GAUSSIAN':
                GAUinp('cgau-s-hf.com',Q,'supHFcca')
            else:
                out_error('ERROR: No program to run Fragpy - set param.program in opt_param.py\n') 

    # ---> Bonded Fragment
    for B1 in frags.bon1:
        if param.program.upper() == 'PSI4':
            wpsi_b(B1,False)
        elif param.program.upper() == 'MOLPRO':
            wmol_b(B1,False)
        elif param.program.upper() == 'GAUSSIAN':
            wgau_b(B1,False)
        else:
            out_error('ERROR: No program to run Fragpy - set param.program in opt_param.py\n') 
    if param.coarse:
        for B1 in frags.bon2:            
            if param.program.upper() == 'PSI4':
                wpsi_b(B1,True)
            elif param.program.upper() == 'MOLPRO':
                wmol_b(B1,True)
            elif param.program.upper() == 'GAUSSIAN':
                wgau_b(B1,True) 
            else:
                out_error('ERROR: No program to run Fragpy - set param.program in opt_param.py\n') 
        for B1 in frags.nbon2:
            if param.program.upper() == 'PSI4':
                wpsi_nb_mb(B1,True)
            elif param.program.upper() == 'MOLPRO':
                wmol_nb_mb(B1,True)
            elif param.program.upper() == 'GAUSSIAN':
                wgau_nb_mb(B1,True)
            else:
                out_error('ERROR: No program to run Fragpy - set param.program in opt_param.py\n') 
            
    # ---> Non-bonded Fragment
    if param.con == 'fnb':
        if param.l1o1:
            for B1 in frags.nbon1:
                if param.program.upper() == 'PSI4':
                    wpsi_nb_mb(B1,False)
                elif param.program.upper() == 'MOLPRO':
                    wmol_nb_mb(B1,False)
                elif param.program.upper() == 'GAUSSIAN':
                    wgau_nb_mb(B1,False)
                else:
                    out_error('ERROR: No program to run Fragpy - set param.program in opt_param.py\n') 
                
    inpf = INPsort(Q,frags)
    return(inpf,frags)

def INPsort(Q,frags):
    
    #--- Weights for Input sorting
    swei = 6.0   # HF sup
    bwei = 4.0   # Bonded
    dwei = 2.0   # Dimer Nbon

    # Initial Input order
    # sup-Csup-Bon-Cbon-nbon-Cnbon
    inpf = []
    elist = []
    tmpe = 0
    if not param.coarse:
        if param.program.upper() == 'PSI4':
            inpf.append('psi-s-hf.in')
        elif param.program.upper() == 'MOLPRO':
            inpf.append('mol-s-hf.in')
        elif param.program.upper() == 'GAUSSIAN':
            inpf.append('gau-s-hf.com')
        for i in Q:
            tmpe += Nelectron[i.atm] * swei
        elist.append(tmpe); tmpe = 0
    if param.coarse and param.coarse_level == 2:
        if param.program.upper() == 'PSI4':
            inpf.append('cpsi-s-hf.in')
        elif param.program.upper() == 'MOLPRO':
            inpf.append('cmol-s-hf.in')
        elif param.program.upper() == 'GAUSSIAN':
            inpf.append('cgau-s-hf.com')
        for i in Q:
            tmpe += Nelectron[i.atm] * swei
        elist.append(tmpe); tmpe = 0
    for i in frags.bon1:
        if not param.program.upper() == 'GAUSSIAN':
            inpf.append(i.name+'.in')
        else:
            inpf.append(i.name+'.com')
        for j in i.coord:
            tmpe += Nelectron[j.atm] * bwei
        elist.append(tmpe); tmpe = 0
    if param.coarse:
        for i in frags.bon2:
            if not param.program.upper() == 'GAUSSIAN':
                inpf.append(i.name+'.in')
            else:
                inpf.append(i.name+'.com')
            for j in i.coord:
                tmpe += Nelectron[j.atm] * bwei
            elist.append(tmpe); tmpe = 0
        if param.l1o1:
            for i in frags.nbon2:
                dnam = 'c'+str(i.frag[0])+str(i.frag[1])
                dele = 0
                if not param.program.upper() == 'GAUSSIAN':
                    inpf.append(dnam+'-'+str(i.frag[0])+'.in')
                else:
                    inpf.append(dnam+'-'+str(i.frag[0])+'.com')
                for j in i.coord.A:
                    tmpe += Nelectron[j.atm]; dele += Nelectron[j.atm]*dwei
                elist.append(tmpe); tmpe = 0
                if not param.program.upper() == 'GAUSSIAN':
                    inpf.append(dnam+'-'+str(i.frag[1])+'.in')
                else:
                    inpf.append(dnam+'-'+str(i.frag[1])+'.com')
                for j in i.coord.B:
                    tmpe += Nelectron[j.atm]; dele += Nelectron[j.atm]*dwei
                elist.append(tmpe); tmpe = 0
                if not param.program.upper() == 'GAUSSIAN':
                    inpf.append(dnam+'.in'); elist.append(dele)
                else:
                    inpf.append(dnam+'.com'); elist.append(dele)
    if param.con == 'fnb':
        if param.l1o1:
            for i in frags.nbon1:
                dnam = str(i.frag[0])+str(i.frag[1])
                dele = 0
                if not param.program.upper() == 'GAUSSIAN':
                    inpf.append(dnam+'-'+str(i.frag[0])+'.in')
                else:
                    inpf.append(dnam+'-'+str(i.frag[0])+'.com')
                for j in i.coord.A:
                    tmpe += Nelectron[j.atm]; dele += Nelectron[j.atm]*dwei
                elist.append(tmpe); tmpe = 0
                if not param.program.upper() == 'GAUSSIAN':
                    inpf.append(dnam+'-'+str(i.frag[1])+'.in')
                else:
                    inpf.append(dnam+'-'+str(i.frag[1])+'.com')
                for j in i.coord.B:
                    tmpe += Nelectron[j.atm]; dele += Nelectron[j.atm]*dwei
                elist.append(tmpe); tmpe = 0
                if not param.program.upper() == 'GAUSSIAN':
                    inpf.append(dnam+'.in'); elist.append(dele)
                else:
                    inpf.append(dnam+'.com'); elist.append(dele)
    elist, inpf = (list(t) for t in zip(*sorted(zip(elist,inpf),
                                        reverse=True)))

    return inpf

def PSIinp(f,L1,check):
    pw = open(f,'w')
    psihead(pw)
    for i in range(len(L1)):
        L2 = L1[i]
        pw.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(L2.atm,L2.qx*Ang,L2.qy*Ang,L2.qz*Ang))
    psibot1(pw)
    if check == 'supHF':
        psibot2(pw,param.basis,param.basis_jk,param.basis_mp2)
        psibotscf(pw)
    elif check == 'supHFcca':
        psibot2(pw,param.cbasis,param.cbasis_jk,param.cbasis_mp2)
        psibotscfa(pw)
    elif check == 'supMP2':
        psibot2(pw,param.basis,param.basis_jk,param.basis_mp2)
        psibotmp2(pw)        
    pw.close()
    
def wpsi_b(Bon,chk):
    if not chk:
        pw = open(Bon.name+'.in','w')
    else:
        pw = open(Bon.name+'.in','w')
    psihead(pw)
    for i in Bon.coord:
        pw.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
    for i in Bon.cap:
        pw.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
    for i in Bon.dum:
        pw.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format('@'+i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
    psibot1(pw)
    if chk and param.coarse_level == 2:
        psibot2(pw,param.cbasis,param.cbasis_jk,param.cbasis_mp2)
        psibotscfa(pw)
    
    psibot2(pw,param.basis,param.basis_jk,param.basis_mp2)
    psibotscf(pw)
    if not chk:
        psibotmp2(pw)
    pw.close()

def wpsi_nb_mb(Nbon,chk):
    name = str(Nbon.frag[0])+str(Nbon.frag[1])
    if not chk:
        w1 = open(name+'.in','w')
        w2 = open(name+'-'+str(Nbon.frag[0])+'.in','w')
        w3 = open(name+'-'+str(Nbon.frag[1])+'.in','w')
    else:
        w1 = open('c'+name+'.in','w')
        w2 = open('c'+name+'-'+str(Nbon.frag[0])+'.in','w')
        w3 = open('c'+name+'-'+str(Nbon.frag[1])+'.in','w')

    psihead(w1); psihead(w2); psihead(w3)
    
    for i in Nbon.coord.A:
        w1.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
        w2.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))

    for i in Nbon.cap.A:
        w2.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))

    for i in Nbon.coord.B:
        w1.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
        w2.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format('@'+i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
        w3.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))

    for i in Nbon.cap.A:
        w1.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))

    for i in Nbon.cap.B:
        w1.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
        w2.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format('@'+i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
        w3.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))

    for i in Nbon.coord.A:
        w3.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format('@'+i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
    psibot1(w1); psibot1(w2); psibot1(w3)
    if chk and param.coarse_level == 2:
        psibot2(w1,param.cbasis,param.cbasis_jk,param.cbasis_mp2)
        psibot2(w2,param.cbasis,param.cbasis_jk,param.cbasis_mp2)
        psibot2(w3,param.cbasis,param.cbasis_jk,param.cbasis_mp2)
        psibotscfa(w1); psibotscfa(w2); psibotscfa(w3)

    psibot2(w1,param.basis,param.basis_jk,param.basis_mp2)
    psibot2(w2,param.basis,param.basis_jk,param.basis_mp2) 
    psibot2(w3,param.basis,param.basis_jk,param.basis_mp2)
    psibotscf(w1); psibotscf(w2); psibotscf(w3)
    if not chk:
        psibotmp2(w1); psibotmp2(w2); psibotmp2(w3)
    w1.close(); w2.close(); w3.close()
    
def psihead(psir):
    psir.write('memory '+str(param.memory)+' mb\n')
    psir.write('\n')
    psir.write('molecule cp {\n')
    psir.write('0 1\n')

def psibot1(psir):
    psir.write('symmetry c1\n')
    psir.write('}\n')
    psir.write('\n')

def psibot2(psir,bas1,bas2,bas3):
    psir.write('\n')
    psir.write('set {\n')
    psir.write('   basis '+bas1+'\n')
    psir.write('   df_basis_scf '+bas2+'\n')
    psir.write('   df_basis_mp2 '+bas3+'\n')
    psir.write('   scf_type df\n')
    psir.write('   guess sad\n')
    psir.write('   freeze_core true\n')
    psir.write('   e_convergence      '+param.E_conv+'\n')
    psir.write('   d_convergence      '+param.D_conv+'\n')
    psir.write('}\n')
    psir.write('\n')
    psir.write('clean() \n')

def psibotscf(psir):
    psir.write('G = gradient(\''+param.cmethod+'\')\n')      
    psir.write('import numpy\n')
    psir.write('G = numpy.array(G)\n')
    psir.write('print_out(\'START SCF GRADIENT\\n\')\n')
    psir.write('for i in G:\n')
    psir.write('    print_out(\'{:<.14e}      {:<.14e}   \
    {:<.14e}\\n \'.format(i[0],i[1],i[2]))\n')
    psir.write('print_out(\'STOP SCF GRADIENT\\n\')\n')
    psir.write('\n')
    psir.write('clean() \n')

def psibotscfa(psir):
    psir.write('G = gradient(\''+param.cmethod+'\')\n')    
    psir.write('import numpy\n')
    psir.write('G = numpy.array(G)\n')
    psir.write('print_out(\'START STO-3G SCF GRADIENT\\n\')\n')
    psir.write('for i in G:\n')
    psir.write('    print_out(\'{:<.14e}      {:<.14e}   \
    {:<.14e}\\n \'.format(i[0],i[1],i[2]))\n')
    psir.write('print_out(\'STOP STO-3G SCF GRADIENT\\n\')\n')
    psir.write('\n')
    psir.write('clean() \n')


def psibotmp2(psir):
    psir.write('G1 = gradient(\''+param.method+'\')\n')
    psir.write('import numpy\n')
    psir.write('G1 = numpy.array(G1)\n')
    psir.write('print_out(\'START MP2 GRADIENT\\n\')\n')
    psir.write('for i in G1:\n')
    psir.write('    print_out(\'{:<.14e}      {:<.14e}     \
    {:<.14e}\\n \'.format(i[0],i[1],i[2]))\n')
    psir.write('print_out(\'STOP MP2 GRADIENT\\n\')\n')
    psir.write('\n')
    psir.write('clean() \n')


# ---> Molpro input files
def MOLinp(f,L1,check):
    pw = open(f,'w')
    molhead(pw)
    for i in range(len(L1)):
        L2 = L1[i]
        pw.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(L2.atm,L2.qx*Ang,L2.qy*Ang,L2.qz*Ang))
    pw.write('}\n')
    if check == 'supHF':
        molbasis(pw,param.basis,param.basis_jk,param.basis_mp2)
        molc(pw)
    elif check == 'supHFcca':
        molbasis(pw,param.cbasis,param.cbasis_jk,param.cbasis_mp2)
        molc(pw)
    elif check == 'supMP2':
        molbasis(pw,param.basis,param.basis_jk,param.basis_mp2)
        molc(pw)
        molf(pw)
    pw.close()

def wmol_b(Bon,chk):
    if not chk:
        pw = open(Bon.name+'.in','w')
    else:
        pw = open(Bon.name+'.in','w')

    molhead(pw)
    for i in Bon.coord:
        pw.write('{:4s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm+str(i.n),i.qx*Ang,i.qy*Ang,i.qz*Ang))
    for i in Bon.cap:
        pw.write('{:4s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm+str(i.n),i.qx*Ang,i.qy*Ang,i.qz*Ang))
    for i in Bon.dum:
        pw.write('{:4s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm+str(i.n),i.qx*Ang,i.qy*Ang,i.qz*Ang))
    pw.write('}\n')
    if chk and param.coarse_level == 2:
        molbasis(pw,param.cbasis,param.cbasis_jk,param.cbasis_mp2)
        pw.write('text,START COARSE FRAGMENT\n')
        moldum(pw,Bon.dum)
        molc(pw)
        pw.write('text,STOP COARSE FRAGMENT\n')

    molbasis(pw,param.basis,param.basis_jk,param.basis_mp2)
    pw.write('text,START DENSE FRAGMENT\n')
    moldum(pw,Bon.dum)
    molc(pw)
    pw.write('text,STOP DENSE FRAGMENT\n')
    if not chk:
        pw.write('text,START MAIN FRAGMENT\n')
        molf(pw)
        pw.write('text,STOP MAIN FRAGMENT\n')
    pw.close()

def wmol_nb_mb(Nbon,chk):
    name = str(Nbon.frag[0])+str(Nbon.frag[1])
    if not chk:
        w1 = open(name+'.in','w')
        w2 = open(name+'-'+str(Nbon.frag[0])+'.in','w')
        w3 = open(name+'-'+str(Nbon.frag[1])+'.in','w')
    else:
        w1 = open('c'+name+'.in','w')
        w2 = open('c'+name+'-'+str(Nbon.frag[0])+'.in','w')
        w3 = open('c'+name+'-'+str(Nbon.frag[1])+'.in','w')

    molhead(w1); molhead(w2); molhead(w3)
    
    for i in Nbon.coord.A:
        w1.write('{:4s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm+str(i.n),i.qx*Ang,i.qy*Ang,i.qz*Ang))
        w2.write('{:4s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm+str(i.n),i.qx*Ang,i.qy*Ang,i.qz*Ang))

    for i in Nbon.cap.A:
        w2.write('{:4s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm+str(i.n)+'A',i.qx*Ang,i.qy*Ang,i.qz*Ang))

    for i in Nbon.coord.B:
        w1.write('{:4s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm+str(i.n),i.qx*Ang,i.qy*Ang,i.qz*Ang))
        w2.write('{:4s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm+str(i.n),i.qx*Ang,i.qy*Ang,i.qz*Ang))
        w3.write('{:4s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm+str(i.n),i.qx*Ang,i.qy*Ang,i.qz*Ang))

    for i in Nbon.cap.A:
        w1.write('{:4s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm+str(i.n),i.qx*Ang,i.qy*Ang,i.qz*Ang))

    for i in Nbon.cap.B:
        w1.write('{:4s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm+str(i.n),i.qx*Ang,i.qy*Ang,i.qz*Ang))
        w2.write('{:4s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm+str(i.n)+'B',i.qx*Ang,i.qy*Ang,i.qz*Ang))
        w3.write('{:4s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm+str(i.n)+'B',i.qx*Ang,i.qy*Ang,i.qz*Ang))

    for i in Nbon.coord.A:
        w3.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm+str(i.n),i.qx*Ang,i.qy*Ang,i.qz*Ang))
    w1.write('}\n'); w2.write('}\n'); w3.write('}\n')

    if chk and param.coarse_level == 2:
        molbasis(w1,param.cbasis,param.cbasis_jk,param.cbasis_mp2)
        molbasis(w2,param.cbasis,param.cbasis_jk,param.cbasis_mp2)
        molbasis(w3,param.cbasis,param.cbasis_jk,param.cbasis_mp2)
        w1.write('text,START COARSE FRAGMENT\n')
        w1.write('text,START COARSE FRAGMENT\n')
        w1.write('text,START COARSE FRAGMENT\n')
        moldum1(w2,Nbon.coord.B+Nbon.cap.B,'B'); moldum1(w3,Nbon.coord.A+Nbon.cap.A,'A')
        molc(w1); molc(w2); molc(w3)
        w1.write('text,STOP COARSE FRAGMENT\n')
        w1.write('text,STOP COARSE FRAGMENT\n')
        w1.write('text,STOP COARSE FRAGMENT\n')

    molbasis(w1,param.basis,param.basis_jk,param.basis_mp2)
    molbasis(w2,param.basis,param.basis_jk,param.basis_mp2) 
    molbasis(w3,param.basis,param.basis_jk,param.basis_mp2)
    w1.write('text,START DENSE FRAGMENT\n')
    w1.write('text,START DENSE FRAGMENT\n')
    w1.write('text,START DENSE FRAGMENT\n')
    moldum1(w2,Nbon.coord.B+Nbon.cap.B,'B'); moldum1(w3,Nbon.coord.A+Nbon.cap.A,'A')
    molc(w1); molc(w2); molc(w3)
    w1.write('text,STOP DENSE FRAGMENT\n')
    w1.write('text,STOP DENSE FRAGMENT\n')
    w1.write('text,STOP DENSE FRAGMENT\n')
    if not chk:
        w1.write('text,START MAIN FRAGMENT\n')
        w1.write('text,START MAIN FRAGMENT\n')
        w1.write('text,START MAIN FRAGMENT\n')
        molf(w1); molf(w2); molf(w3)
        w1.write('text,STOP MAIN FRAGMENT\n')
        w1.write('text,STOP MAIN FRAGMENT\n')
        w1.write('text,STOP MAIN FRAGMENT\n')
    w1.close(); w2.close(); w3.close()

def molhead(f1):
    f1.write('gdirect\n')
    f1.write('gthresh,energy='+param.E_conv.replace('e','d')+',orbital='+param.O_conv.replace('e','d')+'\n')
    f1.write('memory,'+str(param.memory)+',m\n')
    f1.write('\n')
    f1.write('angstrom\n')
    f1.write('geometry={\n')

def moldum(f1,d1):
    f1.write('\n')
    f1.write('dummy')
    for i in d1:
        f1.write(','+i.atm+str(i.n))
    f1.write('\n')
    f1.write('\n')

def moldum1(f1,d1,d2):
    f1.write('\n')
    f1.write('dummy')
    for i in d1:
        if i.n:
            f1.write(','+i.atm+str(i.n))
        else:
            f1.write(','+i.atm+str(i.n)+d2)
    f1.write('\n')
    f1.write('\n')

def molbasis(f1,bas,jbas,mbas):
    f1.write('\n')
    f1.write('basis={\n')
    f1.write('set,orbital\n')
    f1.write('default,'+bas+'\n')
    f1.write('set,jkfit\n')
    f1.write('default,'+jbas+'/jkfit\n')
    f1.write('set,mp2fit\n')
    f1.write('default,'+mbas+'/mp2fit\n')
    f1.write('}')
    f1.write('\n')

def molc(f1):
    f1.write('\n')
    f1.write('{df-'+param.cmethod+',basis=jkfit}\n')
    f1.write('forces\n')
    f1.write('\n')

def molf(f1):
    f1.write('\n')
    f1.write('{df-'+param.method+',basis=mp2fit}\n')
    f1.write('forces\n')
    f1.write('\n')


# ---> Gaussian inputs
def GAUinp(f,L1,check):
    pw = open(f,'w')
    
    if check == 'supHF':
        gauhead(pw,param.cmethod,param.basis,'SUP Correct\n')
    elif check == 'supHFcca':
        gauhead(pw,param.cmethod,param.cbasis,'SUP Correct\n')
    elif check == 'supMP2':
        gauhead(pw,param.method,param.basis,'SUP Fragment\n')
    for i in range(len(L1)):
        L2 = L1[i]
        pw.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(L2.atm,L2.qx*Ang,L2.qy*Ang,L2.qz*Ang))
    pw.write('\n')
    pw.write('\n')       
    pw.close()
    
def wgau_b(Bon,chk):
    if not chk:
        pw = open(Bon.name+'.com','w')
    else:
        pw = open(Bon.name+'.com','w')

    if chk and param.coarse_level == 2:
        gauhead(pw,param.cmethod,param.cbasis,'COARSE FRAGMENTATION\n')
        for i in Bon.coord:
            pw.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                     format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
        for i in Bon.cap:
            pw.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                     format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
        for i in Bon.dum:
            pw.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                     format(i.atm+'-Bq',i.qx*Ang,i.qy*Ang,i.qz*Ang))
        pw.write('\n')
        pw.write('--LINK1--\n')

    gauhead(pw,param.cmethod,param.basis,'DENSE FRAGMENTATION\n')
    for i in Bon.coord:
        pw.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
    for i in Bon.cap:
        pw.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
    for i in Bon.dum:
        pw.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm+'-Bq',i.qx*Ang,i.qy*Ang,i.qz*Ang))
    pw.write('\n'); pw.write('\n')
    
    if not chk:
        pw.write('--LINK1--\n')
        gauhead(pw,param.method,param.basis,'MAIN FRAGMENTATION\n')
        for i in Bon.coord:
            pw.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                     format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
        for i in Bon.cap:
            pw.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                     format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
        for i in Bon.dum:
            pw.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                     format(i.atm+'-Bq',i.qx*Ang,i.qy*Ang,i.qz*Ang))
        pw.write('\n'); pw.write('\n')
    pw.close()

def wgau_nb_mb(Nbon,chk):
    name = str(Nbon.frag[0])+str(Nbon.frag[1])
    if not chk:
        w1 = open(name+'.com','w')
        w2 = open(name+'-'+str(Nbon.frag[0])+'.com','w')
        w3 = open(name+'-'+str(Nbon.frag[1])+'.com','w')
    else:
        w1 = open('c'+name+'.com','w')
        w2 = open('c'+name+'-'+str(Nbon.frag[0])+'.com','w')
        w3 = open('c'+name+'-'+str(Nbon.frag[1])+'.com','w')

    if chk and param.coarse_level == 2:
        gauhead(w1,param.cmethod,param.cbasis,'COARSE FRAGMENTATION\n')
        gauhead(w2,param.cmethod,param.cbasis,'COARSE FRAGMENTATION\n')
        gauhead(w3,param.cmethod,param.cbasis,'COARSE FRAGMENTATION\n')
        gau_nb_geom(w1,w2,w3,Nbon)
        w1.write('\n'); w2.write('\n'); w3.write('\n')
        w1.write('--LINK1--\n'); w2.write('--LINK1--\n'); w3.write('--LINK1--\n')

    gauhead(w1,param.cmethod,param.basis,'DENSE FRAGMENTATION\n')
    gauhead(w2,param.cmethod,param.basis,'DENSE FRAGMENTATION\n')
    gauhead(w3,param.cmethod,param.basis,'DENSE FRAGMENTATION\n')
    gau_nb_geom(w1,w2,w3,Nbon)
    w1.write('\n'); w2.write('\n'); w3.write('\n')
    w1.write('\n'); w2.write('\n'); w3.write('\n')
    if not chk:
        w1.write('--LINK1--\n'); w2.write('--LINK1--\n'); w3.write('--LINK1--\n')
        gauhead(w1,param.method,param.basis,'MAIN FRAGMENTATION\n')
        gauhead(w2,param.method,param.basis,'MAIN FRAGMENTATION\n')
        gauhead(w3,param.method,param.basis,'MAIN FRAGMENTATION\n')
        gau_nb_geom(w1,w2,w3,Nbon)
        w1.write('\n'); w2.write('\n'); w3.write('\n')
        w1.write('\n'); w2.write('\n'); w3.write('\n')
    w1.close(); w2.close(); w3.close()

def gau_nb_geom(w1,w2,w3,Nbon):

    for i in Nbon.coord.A:
        w1.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
        w2.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))

    for i in Nbon.cap.A:
        w2.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))

    for i in Nbon.coord.B:
        w1.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
        w2.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm+'-Bq',i.qx*Ang,i.qy*Ang,i.qz*Ang))
        w3.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))

    for i in Nbon.cap.A:
        w1.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))

    for i in Nbon.cap.B:
        w1.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))
        w2.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm+'-Bq',i.qx*Ang,i.qy*Ang,i.qz*Ang))
        w3.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm,i.qx*Ang,i.qy*Ang,i.qz*Ang))

    for i in Nbon.coord.A:
        w3.write('{:3s}   {:15.10f}   {:15.10f}   {:15.10f}\n'.\
                 format(i.atm+'-Bq',i.qx*Ang,i.qy*Ang,i.qz*Ang))
    
def gauhead(psir,method,basis,comment):
    psir.write('%nprocshared='+str(param.procs)+'\n')
    psir.write('%mem='+str(param.memory)+'MB'+'\n')
    psir.write('%NoSave\n')
    psir.write('# force '+method+'/'+basis)
    for i in param.gau_command:
        psir.write(' '+i+' ')
    psir.write('\n')
    psir.write('\n')
    psir.write(comment)
    psir.write('\n')
    psir.write('0 1\n')

class Fragments:

    def __init__(self):
        self.bon1 = [] #bonded
        for i in range(len(param.bondedlist)):
            self.bon1.append(Fprop(param.bondedlist[i],param.bcoef[i],True))
        if param.coarse:
            self.bon2 = []
            for i in range(len(param.cbondedlist)):
                self.bon2.append(Fprop(param.cbondedlist[i],param.cbcoef[i],False))
        
        self.nbon1 = [] #nbonded
        for i in range(len(param.nbondedlist)):
            self.nbon1.append(Fnprop(param.nbondedlist[i],1.0))
        if param.coarse:
            self.nbon2 = []
            for i in range(len(param.cnbondedlist)):
                self.nbon2.append(Fnprop(param.cnbondedlist[i],1.0))
                        
                
    def addcoord(self,Q):
        cout = 0
        for i1,i2 in param.bondedlist:
            empty = []
            for j in Q:
                if j.grp == i1 or j.grp == i2 :
                    empty.append(j)
            self.bon1[cout].addcoord(empty)
            cout += 1
        if param.coarse:
            cout = 0
            for i1,i2 in param.cbondedlist:
                empty = []
                for j in Q:
                    if j.cgrp == i1 or j.cgrp == i2:
                        empty.append(j)
                self.bon2[cout].addcoord(empty)
                cout += 1

    def addcoord_nb(self,Q):
        cout = 0
        for i1,i2 in param.nbondedlist:
            tmp1 = []; tmp2 = []
            for j in Q:
                if j.grp == i1:
                    tmp1.append(j)
                elif j.grp == i2:
                    tmp2.append(j)
            self.nbon1[cout].coord.addAB(tmp1,tmp2)
            cout += 1
        if param.coarse:
            cout = 0
            for i1,i2 in param.cnbondedlist:
                tmp1 = []; tmp2 = []
                for j in Q:
                    if j.cgrp == i1:
                        tmp1.append(j)
                    elif j.cgrp == i2:
                        tmp2.append(j)
                self.nbon2[cout].coord.addAB(tmp1,tmp2)
                cout += 1

    def addcaps(self,Q):
        cout = 0
        for i1,i2 in param.bondedlist:
            tmp = []; tmp1 = []
            if 0 < i1-1:
                c1,c2 = HCAP(self.bon1[cout].coord,Q,i1-1,1.0)
                tmp.append(c1); tmp1.append(c2)            
            if param.nfrag+1 > i2+1:
                if i2: # overlapped fragment
                    c1,c2 = HCAP(self.bon1[cout].coord,Q,i2+1,1.0)
                    tmp.append(c1); tmp1.append(c2)
                else: # overlapping region
                    c1,c2 = HCAP(self.bon1[cout].coord,Q,i1+1,1.0)
                    tmp.append(c1); tmp1.append(c2)
            self.bon1[cout].addcaps(tmp)
            self.bon1[cout].adddchk(tmp1)
            cout += 1
            
        if param.coarse:
            cout = 0
            for i1,i2 in param.cbondedlist:
                tmp = []; tmp1 = []
                if 0 < i1-1:
                    c1,c2 = HCAPcc(self.bon2[cout].coord,Q,i1-1,1.0)
                    tmp.append(c1); tmp1.append(c2)                
                if param.cnfrag+1 > i2+1:
                    if i2: # overlapped fragment
                        c1,c2 = HCAPcc(self.bon2[cout].coord,Q,i2+1,1.0)
                        tmp.append(c1); tmp1.append(c2)
                    else: # overlapping region
                        c1,c2 = HCAPcc(self.bon2[cout].coord,Q,i1+1,1.0)
                        tmp.append(c1); tmp1.append(c2)
                self.bon2[cout].addcaps(tmp)
                self.bon2[cout].adddchk(tmp1)
                cout += 1
                
    def addcaps_nb(self,Q):
        cout = 0
        nbscale = param.nbscale
        for i1,i2 in param.nbondedlist:
            tmp1 = []; tmp2 = []
            if 0 < i1-1:
                c1,c2 = HCAP(self.nbon1[cout].coord.A,Q,i1-1,nbscale)
                tmp1.append(c1)
            if param.nfrag+1 > i1+1:
                c1,c2 = HCAP(self.nbon1[cout].coord.A,Q,i1+1,nbscale)
                tmp1.append(c1)
            if 0 < i2-1:
                c1,c2 = HCAP(self.nbon1[cout].coord.B,Q,i2-1,nbscale)
                tmp2.append(c1)
            if param.nfrag+1 > i2+1:
                c1,c2 = HCAP(self.nbon1[cout].coord.B,Q,i2+1,nbscale)
                tmp2.append(c1)
            self.nbon1[cout].cap.addAB(tmp1,tmp2)
            cout += 1

        if param.coarse:
            cout = 0
            for i1,i2 in param.cnbondedlist:
                tmp1 = []; tmp2 = []
                if 0 < i1-1:
                    c1,c2 = HCAPcc(self.nbon2[cout].coord.A,Q,i1-1,nbscale)
                    tmp1.append(c1)
                if param.cnfrag+1 > i1+1:
                    c1,c2 = HCAPcc(self.nbon2[cout].coord.A,Q,i1+1,nbscale)
                    tmp1.append(c1)
                if 0 < i2-1:
                    c1,c2 = HCAPcc(self.nbon2[cout].coord.B,Q,i2-1,nbscale)
                    tmp2.append(c1)
                if param.cnfrag+1 > i2+1:
                    c1,c2 = HCAPcc(self.nbon2[cout].coord.B,Q,i2+1,nbscale)
                    tmp2.append(c1)
                self.nbon2[cout].cap.addAB(tmp1,tmp2)
                cout += 1
            
            

    def adddums(self,Q):
        if not param.dumfrag:
            cout = 0
            for i1,i2 in param.bondedlist:
                tmp1 = []
                for i in Q:
                    if (not i in self.bon1[cout].dchk and
                        not i.grp == i1 and not i.grp == i2):
                        for j in self.bon1[cout].coord:
                            d = norm(subtract(i.q,j.q).tolist())
                            if d <= param.bonthres/Ang:
                                if not i in tmp1:
                                    tmp1.append(i)
                self.bon1[cout].adddums(tmp1)
                cout += 1
            if param.coarse:
                cout = 0
                for i1,i2 in param.cbondedlist:
                    tmp1 = []
                    for i in Q:
                        if (not i in self.bon2[cout].dchk and
                            not i.cgrp == i1 and not i.cgrp == i2):
                            for j in self.bon2[cout].coord:
                                d = norm(subtract(i.q,j.q).tolist())
                                if d <= param.cbonthres/Ang:
                                    if not i in tmp1:
                                        tmp1.append(i)
                    self.bon2[cout].adddums(tmp1)
                    cout += 1
        else:
            cout = 0
            tmp2 = param.dumfrag_level
            for i1,i2 in param.bondedlist:
                tmp1 = []
                if 0 < i1-1:
                    Ngrp = []
                    for i in range(i1-tmp2 if i1-tmp2>0 else 1,i1):
                        Ngrp.append(i)
                    tmp1.extend(DumFrag(self.bon1[cout],Q,Ngrp))
                if param.nfrag+1 > i2+1:
                    if i2:
                        Ngrp = []
                        for i in range(i2+1,i2+tmp2+1 if i2+tmp2 <= param.nfrag else param.nfrag+1):
                            Ngrp.append(i)
                        tmp1.extend(DumFrag(self.bon1[cout],Q,Ngrp))
                    else:
                        Ngrp = []
                        for i in range(i1+1,i1+tmp2+1 if i1+tmp2 <= param.nfrag else param.nfrag+1):
                            Ngrp.append(i)
                        tmp1.extend(DumFrag(self.bon1[cout],Q,Ngrp))
                self.bon1[cout].addums(tmp1)
                cout += 1

            if param.coarse:
                cout = 0
                for i1,i2 in param.cbondedlist:
                    tmp1 = []
                    if 0 < i1-1:
                        Ngrp = []
                        for i in range(i1-tmp2 if i1-tmp2>0 else 1,i1):
                            Ngrp.append(i)
                        tmp1.extend(DumFrag1(self.bon2[cout],Q,Ngrp))
                    if param.cnfrag+1 > i2+1:
                        if i2:
                            Ngrp = []
                            for i in range(i2+1,i2+tmp2+1 if i2+tmp2 <= param.cnfrag else param.cnfrag+1):
                                Ngrp.append(i)
                            tmp1.extend(DumFrag1(self.bon2[cout],Q,Ngrp))
                        else:
                            Ngrp = []
                            for i in range(i1+1,i1+tmp2+1 if i1+tmp2 <= param.cnfrag else param.cnfrag+1):
                                Ngrp.append(i)
                            tmp1.extend(DumFrag1(self.bon2[cout],Q,Ngrp))
                    self.bon2[cout].addums(tmp1)
                    cout += 1

class Fprop:
    
    def __init__(self,nam,coef,chk00):
        if chk00:
            self.name = str(nam[0])+'' if nam[1]==0 else str(nam[0])+str(nam[1])
        else:
            self.name = 'c'+str(nam[0])+'' if nam[1]==0 else 'c'+str(nam[0])+str(nam[1])
        self.coef = coef
        self.coord = []
        self.cap = []
        self.dum = []
        self.dchk = []
        
    def addcoord(self,f):
        self.coord = f
    
    def addcaps(self,f):
        self.cap = f
        
    def adddchk(self,f):
        self.dchk = f

    def adddums(self,f):
        self.dum = f
                    

class Fnprop:
    
    def __init__(self,nam,coef):
        self.frag = nam
        self.coef = coef
        self.coord = Monomer()
        self.cap = Monomer()
                

class Monomer:
    def __init__(self):
        self.A = []
        self.B = []
        
    def addAB(self,a,b):
        self.A = a
        self.B = b

def DumFrag(L1,L2,Ngrp):
    tmp = []
    for i in L2:
        if not i in L1.dchk:
            if i.grp in Ngrp:
                tmp.append(i)
    return tmp
            
            
def DumFrag1(L1,L2,Ngrp):
    tmp = []
    for i in L2:
        if not i in L1.dchk:
            if i.cgrp in Ngrp:
                tmp.append(i)
    return tmp                
                

def HCAP(L1,L2,n2,scale):
    min = 1000.0
    for j in L1:
        for k in L2:
            if k.grp == n2:
                d = norm(subtract(j.q,k.q).tolist())
                if d < min:
                    min = d
                    mini = [j,k]
    q1 = mini[0]; q2 = mini[1]
    tb1 = braggr[q1.atm]+braggr['H']
    tb2 = braggr[q1.atm]+braggr[q2.atm]
    scalb = tb1/tb2*scale
    x1 = q1.qx + ((q2.qx-q1.qx)*scalb)
    x2 = q1.qy + ((q2.qy-q1.qy)*scalb)
    x3 = q1.qz + ((q2.qz-q1.qz)*scalb)
    return(Coords(0,[x1,x2,x3],'H',q1.grp,q1.cgrp),q2)
                
def HCAPcc(L1,L2,n2,scale):
    min = 1000.0
    for j in L1:
        for k in L2:
            if k.cgrp == n2:
                d = norm(subtract(j.q,k.q).tolist())
                if d < min:
                    min = d
                    mini = [j,k]
    q1 = mini[0]; q2 = mini[1]
    tb1 = braggr[q1.atm]+braggr['H']
    tb2 = braggr[q1.atm]+braggr[q2.atm]
    scalb = tb1/tb2*scale
    x1 = q1.qx + ((q2.qx-q1.qx)*scalb)
    x2 = q1.qy + ((q2.qy-q1.qy)*scalb)
    x3 = q1.qz + ((q2.qz-q1.qz)*scalb)
    return(Coords(0,[x1,x2,x3],'H',q1.grp,q1.cgrp),q2)
