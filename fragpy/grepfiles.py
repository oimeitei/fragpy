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
Reads gradients and energies from other codes
PSI4 
Molpro
'''

from fragpy.opt_param import param
import numpy as np


def fraggradpsi(frags):
    if param.con == 's':
        if param.program.upper() == 'PSI4':
            sg,se = smp2_psi('psi-s.out')
        elif param.program.upper() == 'MOLPRO':
            sg,se = smp2_mol('mol-s.out')
        elif param.program.upper() == 'GAUSSIAN':
            sg,se = smp2_gau('gau-s.log')
        else:
            out_error('ERROR: No program to run Fragpy - set param.program in opt_param.py\n')
        empty1 = []
        for i in range(param.natom):
            empty1.extend([sg[i].gx,sg[i].gy,sg[i].gz])

        return(empty1,se) 

    Energy = FRAGMENT('E'); Gradient = FRAGMENT('G')

    if param.DEBUG >= 1:
        out = open('output','a')
    for i in frags.bon1:
        ce,fe,cg,fg = readbon(i,0)
        Gradient.addcB(cg); Gradient.addfB(fg)
        Energy.addfB(i.coef*fe); Energy.addcB(i.coef*ce)

    if param.coarse:
        for i in frags.bon2:
            ce,cg = readbon(i,1)
            Energy.addccB(i.coef*ce); Gradient.addccB(cg)
            
            if param.coarse_level == 2:
                ce,cg = readbon(i,2)
                Energy.addccaB(i.coef*ce); Gradient.addccaB(cg)

        # Read NOTE1 at the end of file!
        if param.l1o1:
            for i in frags.nbon2:
                ce,cg = readnbon(i,1)
                Gradient.addccNB(cg); Energy.addccNB(ce)

                if param.coarse_level == 2:
                    ce,cg = readnbon(i,2)
                    Gradient.addccaNB(cg); Energy.addccaNB(ce)
            
    
    if not param.con=='fb':
        if param.l1o1:
            for i in frags.nbon1:
                ce,fe,cg,fg = readnbon(i,0)
                Gradient.addfNB(fg); Gradient.addcNB(cg)
                Energy.addfNB(fe); Energy.addcNB(ce)

    fgcount=0

    if param.DEBUG >= 1:
        out.write('\n'); out.write(' Correction increments used : \n')

    if not param.coarse:
        if param.program.upper() == 'PSI4':
            hf123glist,hf123e=readgfrag_spsi('psi-s-hf')
        elif param.program.upper() == 'MOLPRO':
            hf123glist,hf123e=readgfrag_spsi('mol-s-hf')
        elif param.program.upper() == 'GAUSSIAN':
            hf123glist,hf123e=readgfrag_spsi('gau-s-hf')
        else:
            out_error('ERROR: No program to run Fragpy - set param.program in opt_param.py\n')
        Energy.sup1 = hf123e; Gradient.sup1 = hf123glist

        if param.DEBUG >= 1:
            out.write(' Sup - Bon({:>3}-grps) - NonBon({:>3}-grps) at HF/{:<5} \n'.
                      format(param.nfrag,param.nfrag,param.basis))

    elif param.coarse:
        Gradient.fragmentationcc()
        Energy.correct1()
        if param.DEBUG >= 1:
            out.write(' (Bon({:>3}-grps) + NonBon({:>3}-grps)) - '\
                      'Bon( {:>3}-grps) - NonBon({:>3}-grps) at HF/{:>5} \n'. \
                      format(param.cnfrag,param.cnfrag,param.nfrag,param.nfrag,param.basis))
            
        if param.coarse_level == 2:
            if param.program.upper() == 'PSI4':
                hf123glista,hf123ea=readgfrag_spsi('cpsi-s-hf')
            elif param.program.upper() == 'MOLPRO':
                hf123glista,hf123ea=readgfrag_spsi('cmol-s-hf')
            elif param.program.upper() == 'GAUSSIAN':
                hf123glista,hf123ea=readgfrag_spsi('cgau-s-hf')
            else:
                out_error('ERROR: No program to run Fragpy - set param.program in opt_param.py\n')
            Energy.sup2 = hf123ea; Gradient.sup2 = hf123glista

            correction2 = Gradient.fragmentationcca()
            Energy.correct2()

            if param.DEBUG >= 1:
                out.write('  + {Sup - Bon({:>3}-grps) - '\
                'NonBon({:>3}-grps)} at HF/{:<5} \n'.
                          format(param.cnfrag,param.cnfrag,param.cbasis))
    
    Gbonf=bonfrag(Gradient.bon.f,0); Gbonc=bonfrag(Gradient.bon.c,0)
    
    if not param.con=='fb':
        Gnbonf=nbonfrag(Gradient.nbon.f,0); Gnbonc=nbonfrag(Gradient.nbon.c,0)
        sglist = fragmentation(Gradient.sup1,Gbonf,Gbonc,Gnbonf,Gnbonc)

        if param.DEBUG >= 1:
            Energy.printfE(out)
        
    elif param.con=='fb':
        if param.DEBUG >= 1:
            out = open('output','a')
            out.write('\n')
            out.write(' Neglecting Nonbonded contribution - Correction terms includes it\n')
            out.write('\n')
        sglist=fragmentation1(Gradient.sup1,Gbonf,Gbonc)

        if param.DEBUG >= 1:
            Energy.printfE(out)

    if param.coarse and param.coarse_level == 2:
        sglist = fragmentationcca(sglist,correction2)

    if param.DEBUG >= 2:
        out.write('\n')
        out.write(' DEBUG set to 2 \n')
        out.write(' Printing gradient informations\n')
        out.write('\n')
        out.write(' Nuclear Gradient Computed with IFM\n')
        out.write('\n')
        out.write(' Atm         X                Y                Z     \n')
        cce = 0
        for i1 in Gbonc: #sglist:
            cce += 1
            out.write('{:2}   {:>15.10f}  {:>15.10f}  {:>15.10f}\n'.
                      format(cce,i1.gx,i1.gy,i1.gz))
        out.write('\n')

    if param.DEBUG >= 1:
        out.write('\n')
        out.close()

    empty1 = []
    for i in range(param.natom):
        empty1.extend([sglist[i].gx,sglist[i].gy,sglist[i].gz])

    return(empty1,Energy.Energy())



def fragmentationc1(file1,file2):
    emptylist= np.zeros((param.natom,3))
    for i in range(param.natom):
        emptylist[i][0]=file1[i][0]+file2[i][0]
        emptylist[i][1]=file1[i][1]+file2[i][1]
        emptylist[i][2]=file1[i][2]+file2[i][2]
    return emptylist

def bonfrag(list1,checkw):
    emptylist=[]
    for i in range(param.natom):
        emptylist.append(Gmethod(i+1,np.zeros(3)))
    for i in range(param.natom):
        bfcount=0
        for j in list1: # list of frag
            for k in j: # list of atoms in frag(j)
                if int(k.n) == i+1:
                    if not checkw:
                        tmp_coef = param.bcoef[bfcount]
                    else:
                        tmp_coef = param.cbcoef[bfcount]

                    emptylist[i].gx += tmp_coef*k.gx
                    emptylist[i].gy += tmp_coef*k.gy
                    emptylist[i].gz += tmp_coef*k.gz
            bfcount+=1
    return emptylist

def nbonfrag(list1,check):
    emptylist=[]
    for i in range(param.natom):
        emptylist.append(Gmethod(i+1,np.zeros(3)))

    for i in range(param.natom):
        nbfcount=0
        for j in list1: #list of frags
            for j1 in range(3):
                for k in j[j1]:
                    if int(k.n) == i+1:
                        if check:
                            t_coef = nbcoef[nbfcount]
                        else:
                            t_coef = 1.0
                        if not j1:
                            t1c = 1.0
                        else:
                            t1c = -1.0
                        emptylist[i].gx += t1c*t_coef*k.gx;
                        emptylist[i].gy += t1c*t_coef*k.gy
                        emptylist[i].gz += t1c*t_coef*k.gz
            nbfcount+=1
    return emptylist


def fragmentation(F1,F2,F3,F4,F5):
    emptylist=[]
    for i in range(param.natom):
        xg = F1[i].gx + (F2[i].gx - F3[i].gx) + (F4[i].gx - F5[i].gx)
        yg = F1[i].gy + (F2[i].gy - F3[i].gy) + (F4[i].gy - F5[i].gy)
        zg = F1[i].gz + (F2[i].gz - F3[i].gz) + (F4[i].gz - F5[i].gz)
        emptylist.append(Gmethod(i+1,[xg,yg,zg]))
    return emptylist

def fragmentationcca(F1,F2):
    emptylist=[]
    for i in range(param.natom):
        xg = F1[i].gx + F2[i].gx; yg = F1[i].gy + F2[i].gy; zg = F1[i].gz + F2[i].gz
        emptylist.append(Gmethod(i+1,xg,yg,zg))
    return emptylist


def fragmentation1(F1,F2,F3):
    emptylist=[]
    for i in range(param.natom):
        xg = F1[i].gx + F2[i].gx - F3[i].gx; yg = F1[i].gy + F2[i].gy - F3[i].gy
        zg = F1[i].gz + F2[i].gz - F3[i].gz
        emptylist.append(Gmethod(i+1,[xg,yg,zg]))
    return emptylist

def readnbon(f1,chk):
    if chk == 0:
        r1 = str(f1.frag[0])+str(f1.frag[1])
        r2 = str(f1.frag[0])+str(f1.frag[1])+'-'+str(f1.frag[0])
        r3 = str(f1.frag[0])+str(f1.frag[1])+'-'+str(f1.frag[1])
    elif chk == 1 or chk == 2:
        r1 = 'c'+str(f1.frag[0])+str(f1.frag[1])
        r2 = 'c'+str(f1.frag[0])+str(f1.frag[1])+'-'+str(f1.frag[0])
        r3 = 'c'+str(f1.frag[0])+str(f1.frag[1])+'-'+str(f1.frag[1])

    Dimer = f1.coord.A + f1.coord.B

    if chk == 0:
        ehf0,emp20,ghf0,gmp20 = readg1o1(r1,Dimer,0)
        ehf1,emp21,ghf1,gmp21 = readg1o1(r2,f1.coord.A,0)
        ehf2,emp22,ghf2,gmp22 = readg1o1(r3,f1.coord.B,0)

    elif chk == 1:
        ehf0,ghf0 = readg1o1(r1,Dimer,1)
        ehf1,ghf1 = readg1o1(r2,f1.coord.A,1)
        ehf2,ghf2 = readg1o1(r3,f1.coord.B,1)

    elif chk == 2:
        ehf0,ghf0 = readg1o1(r1,Dimer,2)
        ehf1,ghf1 = readg1o1(r2,f1.coord.A,2)
        ehf2,ghf2 = readg1o1(r3,f1.coord.B,2)

    if chk == 0:
        emptylist1 = [ghf0,ghf1,ghf2]
        emptylist2 = [gmp20,gmp21,gmp22]
        ehf = ehf0 -ehf1 -ehf2
        emp2 = emp20 -emp21 -emp22
        return(ehf,emp2,emptylist1,emptylist2)

    elif chk == 1 or chk == 2:
        emptylist1 = [ghf0,ghf1,ghf2]
        ehf = ehf0 -ehf1 -ehf2
        return(ehf,emptylist1)

def readg1o1(file1,file2,chck0):
    if param.program.upper() == 'PSI4':
        return nb_mb_psi(open(file1+'.out','r'),\
                         file2,chck0)
    elif param.program.upper() == 'MOLPRO':
        return nb_mb_mol(open(file1+'.out','r'),\
                         file2,chck0)
    elif param.program.upper() == 'GAUSSIAN':
        return nb_mb_gau(open(file1+'.log','r'),\
                         file2,chck0)
    else:
        out_error('ERROR: No program to run Fragpy - set param.program in opt_param.py\n')
        
def readbon(f1,chk): 
    if param.program.upper() == 'PSI4':
        return b_psi(f1,chk)
    elif param.program.upper() == 'MOLPRO':
        return b_mol(f1,chk)
    elif param.program.upper() == 'GAUSSIAN':
        return b_gau(f1,chk)
    else:
        out_error('ERROR: No program to run Fragpy - set param.program in opt_param.py\n')

def readgfrag_spsi(file1):
    if param.program.upper() == 'PSI4':
        return s_psi(file1+'.out')
    elif param.program.upper() == 'MOLPRO':
        return s_mol(file1+'.out')
    elif param.program.upper() == 'GAUSSIAN':
        return s_gau(file1+'.log')
    else:
        out_error('ERROR: No program to run Fragpy - set param.program in opt_param.py\n')

# ----------------------
# Read PSI4 output files

def nb_mb_psi(file1,file2,chck0):
    check1 = 0
    check2 = 0
    check3 = 0
    check4 = 0
    check5 = 0
    count1 = 0
    count2 = 0
    Xcount = 0

    Grad = Gradobj()
    scfbas = 'START SCF GRADIENT'
    if chck0 == 2:
        scfbas = 'START STO-3G SCF GRADIENT'

    for line in file1.readlines():
        if scfbas in line:
            check2 += 1
            continue
        elif 'START MP2 GRADIENT' in line:
            check3 += 1
            continue
        elif 'DF-MP2 Energies' in line:
            check4 = 1
        elif '@DF-RHF Final Energy' in line:
            count2 += 1
            if chck0 == 0:
                if count2 == 2:
                    count2 = 0
                    check5 = 1
            elif chck0 == 1:
                if not param.coarse_level == 2:
                    check5 = 1
                else:
                    if count2 == 2:
                        check5 = 1
            elif chck0 == 2:
                if count2 == 1:
                    check5 = 1
        if check2 == 2:
            if count1+1 > len(file2):
                check2 = 0
                count1 = 0
                continue
            else:
                s = line.split()
                s1 = [float(xx) for xx in s]
                Grad.addc(file2[count1].n,s1[:3])
                count1 += 1

        elif check3 == 2:
            if count1+1 > len(file2):
                check3 = 0
                count1 = 0
                continue
            else:
                s = line.split()
                s1 = [float(xx) for xx in s]
                Grad.addf(file2[count1].n,s1[:3])
                count1 += 1
            
        elif check4 == 1:
            if 'Total Energy' in line:
                s = line.split()
                energymp2p = float(s[-2])
                check4 = 0
        elif check5 == 1:
            if 'Total Energy' in line:
                s = line.split()
                energyscfp = float(s[-1])
                check5 = 0
                
    file1.close()
    if chck0 == 0:
        return(energyscfp,energymp2p,Grad.c,Grad.f)
    elif chck0 == 1 or chck0 == 2:
        return(energyscfp,Grad.c)

def b_psi(f1,chk):
    check1 = 0
    check2 = 0
    check3 = 0
    check4 = 0
    check5 = 0
    count2 = 0
    Xcount = 0


    count1 = 0
    Grad = Gradobj()
    psirr = open(f1.name+'.out','r')
    scfbas = 'START SCF GRADIENT'
    if chk == 2:
        scfbas = 'START STO-3G SCF GRADIENT'

    for line in psirr.readlines():
        
        if scfbas in line:
            check2 += 1
            continue
        elif 'START MP2 GRADIENT' in line:
            check3 += 1
            continue
        elif 'DF-MP2 Energies' in line:
            check4 = 1
        elif '@DF-RHF Final Energy' in line:
            count2 += 1
            if chk == 0:
                if count2 == 2:
                    count2 = 0
                    check5 = 1
            elif chk == 1:
                if not param.coarse_level == 2:
                    check5 = 1
                else:
                    if count2 == 2:
                        check5 = 1
            elif chk == 2:
                if count2 == 1:
                    check5 = 1
        if check2 == 2:
            if count1+1 > len(f1.coord):
                check2 = 0
                count1 = 0
                continue
            else:
                s = line.split()
                s1 = [float(xx) for xx in s]
                Grad.addc(f1.coord[count1].n,s1[:3])
                count1 += 1

        elif check3 == 2:
            if count1+1 > len(f1.coord):
                check3 = 0
                count1 = 0
                continue
            else:
                s = line.split()
                s1 = [float(xx) for xx in s]
                Grad.addf(f1.coord[count1].n,s1[:3])
                count1 += 1
                    
        elif check4 == 1:
            if 'Total Energy' in line:
                s = line.split()
                energymp2p = float(s[-2])
                check4 = 0
        elif check5 == 1:
            if 'Total Energy' in line:
                s = line.split()
                energyscfp = float(s[-1])
                check5 = 0

    psirr.close()
    if chk == 0:
        return(energyscfp,energymp2p,Grad.c,Grad.f)
    elif chk == 1 or chk == 2:
        return(energyscfp,Grad.c)

def s_psi(file1):
    check1 = 0
    check2 = 0
    check3 = 0
    check4 = 0
    check5 = 0
    count1 = 0
    count2 = 0
    Grad = Gradobj()
    Xcount = 0
    psirr = open(file1,'r')
    for line in psirr.readlines():
        if 'START SCF GRADIENT' in line:
            check2 += 1
            continue
        elif '@DF-RHF Final Energy' in line:
            check5 = 1
        if check2 == 2:
            if not 'STOP SCF GRADIENT' in line:
                count1 += 1
                s = line.split()
                s1 = [float(xx) for xx in s]
                Grad.adds(0,s1[:3])
            else:
                check2 = 0
        elif check5 == 1:
            if 'Total Energy' in line:
                s = line.split()
                energyscfp = float(s[-1])
                check5 = 0
    psirr.close()
                
    return(Grad.s,energyscfp)

def smp2_psi(file1):
    check1 = 0
    check2 = 0
    check3 = 0
    check4 = 0
    check5 = 0
    count1 = 0
    count2 = 0
    Grad = []
    Xcount = 0
    psirr = open(file1,'r')
    cout = 1
    for line in psirr.readlines():
        if 'START MP2 GRADIENT' in line:
            check2 += 1
            continue
        elif 'DF-MP2 Energies' in line:
            check5 = 1
        if check2 == 2:
            if not 'STOP MP2 GRADIENT' in line:
                count1 += 1
                s = line.split()
                s1 = [float(xx) for xx in s]
                Grad.append(Gmethod(cout,s1[:3]))
                cout += 1
            else:
                check2 = 0
        elif check5 == 1:
            if 'Total Energy' in line:
                s = line.split()
                energyscfp = float(s[-2])
                check5 = 0
    psirr.close()
                
    return(Grad,energyscfp)

# ---> Molpro

def nb_mb_mol(file1,file2,chck0):
    check1 = 0
    check2 = 0
    check3 = 0
    check4 = 0
    check5 = 0
    count1 = 0
    count2 = 0
    count3 = 0

    Grad = Gradobj()
    scfbas = 'START DENSE FRAGMENT'
    scfbas1 = 'STOP DENSE FRAGMENT'
    if chck0 == 2:
        scfbas = 'START COARSE FRAGMENT'
        scfbas1 = 'STOP COARSE FRAGMENT'

    for line in file1.readlines():
        if scfbas in line:
            check1 += 1; continue
        elif scfbas1 in line:
            check1 -= 1; continue
        elif 'START MAIN FRAGMENT' in line:
            check2 += 1; continue
        elif 'STOP MAIN FRAGMENT' in line:
            check2 -= 1; continue
        elif '!MP2 total energy' in line:
            check3 = 1
        elif '!RHF STATE 1.1 Energy' in line:
            count2 += 1
            if chck0 == 0:
                if count2 == 1:
                    count2 = 0
                    check4 = 1
            elif chck0 == 1:
                if not param.coarse_level == 2:
                    check4 = 1
                else:
                    if count2 == 2:
                        check4 = 1
            elif chck0 == 2:
                if count2 == 1:
                    check4 = 1
            
        if check3 == 1:
            s = line.split()
            energymp2p = float(s[-1])
            check3 = 0
        elif check4 == 1:
            s = line.split()
            energyscfp = float(s[-1])
            check4 = 0

        if check1 == 1:
            if 'Atom          dE/dx               dE/dy               dE/dz' in line:
                count3 = 1
                continue
            elif not count3:
                continue
            if (count3 and line.isspace()):
                continue
            if count1+1 > len(file2):
                check1 = 0
                count1 = 0
                count3 = 0
                continue
            else:
                s = line.split()
                s1 = [float(xx) for xx in s]
                Grad.addc(file2[count1].n,s1[1:])
                count1 += 1

        elif check2 == 1:
            if 'Atom          dE/dx               dE/dy               dE/dz' in line:
                count3 = 1
                continue
            elif not count3:
                continue
            if (count3 and line.isspace()):
                continue
            if count1+1 > len(file2):
                check2 = 0
                count1 = 0
                continue
            else:
                s = line.split()
                s1 = [float(xx) for xx in s]
                Grad.addf(file2[count1].n,s1[1:])
                count1 += 1
                
    file1.close()
    if chck0 == 0:
        return(energyscfp,energymp2p,Grad.c,Grad.f)
    elif chck0 == 1 or chck0 == 2:
        return(energyscfp,Grad.c)


def b_mol(f1,chk):
    check1 = 0
    check2 = 0
    check3 = 0
    check4 = 0
    check5 = 0
    count1 = 0
    count2 = 0
    count3 = 0

    Grad = Gradobj()
    psirr = open(f1.name+'.out','r')

    scfbas = 'START DENSE FRAGMENT'
    scfbas1 = 'STOP DENSE FRAGMENT'
    if chk == 2:
        scfbas = 'START COARSE FRAGMENT'
        scfbas1 = 'STOP COARSE FRAGMENT'

    for line in psirr.readlines():
        
        if scfbas in line:
            check1 += 1; continue
        elif scfbas1 in line:
            check1 -= 1; continue
        elif 'START MAIN FRAGMENT' in line:
            check2 += 1; continue
        elif 'STOP MAIN FRAGMENT' in line:
            check2 -= 1; continue
        elif '!MP2 total energy' in line:
            check3 = 1
        elif '!RHF STATE 1.1 Energy' in line:
            count2 += 1
            if chk == 0:
                if count2 == 1:
                    count2 = 0
                    check4 = 1
            elif chk == 1:
                if not param.coarse_level == 2:
                    check4 = 1
                else:
                    if count2 == 2:
                        check4 = 1
            elif chk == 2:
                if count2 == 1:
                    check4 = 1

        if check3 == 1:
            s = line.split()
            energymp2p = float(s[-1])
            check3 = 0
        elif check4 == 1:
            s = line.split()
            energyscfp = float(s[-1])
            check4 = 0

        if check1 == 1:
            if 'Atom          dE/dx               dE/dy               dE/dz' in line:
                count3 = 1
                continue
            elif not count3:
                continue
            if (count3 and line.isspace()):
                continue
            if count1+1 > len(f1.coord):
                check1 = 0
                count1 = 0
                count3 = 0
                continue
            else:
                s = line.split()
                s1 = [float(xx) for xx in s]
                Grad.addc(f1.coord[count1].n,s1[1:])
                count1 += 1

        elif check2 == 1:
            if 'Atom          dE/dx               dE/dy               dE/dz' in line:
                count3 = 1
                continue
            elif not count3:
                continue
            if (count3 and line.isspace()):
                continue
            if count1+1 > len(f1.coord):
                check2 = 0
                count1 = 0
                continue
            else:
                s = line.split()
                s1 = [float(xx) for xx in s]
                Grad.addf(f1.coord[count1].n,s1[1:])
                count1 += 1

    psirr.close()
    if chk == 0:
        return(energyscfp,energymp2p,Grad.c,Grad.f)
    elif chk == 1 or chk == 2:
        return(energyscfp,Grad.c)

def s_mol(file1):
    check1 = 0
    check2 = 0
    check3 = 0
    check4 = 0
    check5 = 0
    count1 = 0
    count2 = 0
    Grad = Gradobj()
    Xcount = 0
    psirr = open(file1,'r')
    for line in psirr.readlines():
        if 'Atom          dE/dx               dE/dy               dE/dz' in line:
            check1 = 1
            continue
        elif '!RHF STATE 1.1 Energy' in line:
            s = line.split()
            energyscfp = float(s[-1])

        if check1 == 1:
            count1 += 1
            if count1 == 1:
                continue
            if line.isspace():
                check1 = 0
                continue
            s = line.split()
            s1 = [float(xx) for xx in s]
            Grad.adds(0,s1[1:])                
    return(Grad.s,energyscfp)

def smp2_mol(file1):
    check1 = 0
    check2 = 0
    check3 = 0
    check4 = 0
    check5 = 0
    count1 = 0
    count2 = 0
    Grad = []
    Xcount = 0
    cout = 1
    psirr = open(file1,'r')
    for line in psirr.readlines():
        if 'Atom          dE/dx               dE/dy               dE/dz' in line:
            check1 = 1
            continue
        elif '!MP2 total energy' in line:
            s = line.split()
            energyscfp = float(s[-1])

        if check1 == 1:
            count1 += 1
            if count1 == 1:
                continue
            if line.isspace():
                check1 = 0
                continue
            s = line.split()
            s1 = [float(xx) for xx in s]
            Grad.append(Gmethod(cout,s1[1:]))
            cout += 1
                
    return(Grad,energyscfp)

# ---> Gaussian output

def nb_mb_gau(file1,file2,chck0):
    check1 = 0
    check2 = 0
    check3 = 0
    check4 = 0
    check5 = 0
    count1 = 0
    count2 = 0
    count3 = 0

    Grad = Gradobj()
    scfbas = 'DENSE FRAGMENTATION'
    if chck0 == 2:
        scfbas = 'COARSE FRAGMENTATION'

    c1 = 1; c2 = 1
    for line in file1.readlines():
        if scfbas in line:
            check1 += c1*1
            c1 = -1; continue
        elif 'MAIN FRAGMENT' in line:
            check2 += c2*1
            c2 = -1; continue
        elif 'EUMP2' in line:
            check3 = 1
        elif 'SCF Done:  E(RHF)' in line:
            count2 += 1
            if chck0 == 0:
                if count2 == 1:
                    count2 = 0
                    check4 = 1
            elif chck0 == 1:
                if not param.coarse_level == 2:
                    check4 = 1
                else:
                    if count2 == 2:
                        check4 = 1
            elif chck0 == 2:
                if count2 == 1:
                    check4 = 1
            
        if check3 == 1:
            s = line.split()
            energymp2p = float(s[-1].replace('D','e'))
            check3 = 0
        elif check4 == 1:
            s = line.split()
            energyscfp = float(s[4])
            check4 = 0

        if check1 == 1:
            if 'Forces (Hartrees/Bohr)' in line:
                count3 = 1
                continue
            elif not count3:
                continue
            if (count3 and 'Number     Number              X              Y              Z' in line) or\
               (count3 and '------------------' in line):
                continue
            if count1+1 > len(file2):
                check1 = 0
                count1 = 0
                count3 = 0
                continue
            else:
                s = line.split()
                s1 = [-float(xx) for xx in s]
                Grad.addc(file2[count1].n,s1[2:])
                count1 += 1

        elif check2 == 1:
            if 'Forces (Hartrees/Bohr)' in line:
                count3 = 1
                continue
            elif not count3:
                continue
            if (count3 and 'Number     Number              X              Y              Z' in line)or\
               (count3 and '------------------' in line):
                continue
            if count1+1 > len(file2):
                check2 = 0
                count1 = 0
                continue
            else:
                s = line.split()
                s1 = [-float(xx) for xx in s]
                Grad.addf(file2[count1].n,s1[2:])
                count1 += 1
                
    file1.close()
    if chck0 == 0:
        return(energyscfp,energymp2p,Grad.c,Grad.f)
    elif chck0 == 1 or chck0 == 2:
        return(energyscfp,Grad.c)


def b_gau(f1,chk):
    check1 = 0
    check2 = 0
    check3 = 0
    check4 = 0
    check5 = 0
    count1 = 0
    count2 = 0
    count3 = 0

    Grad = Gradobj()
    psirr = open(f1.name+'.log','r')
    scfbas = 'DENSE FRAGMENTATION'
    if chk == 2:
        scfbas = 'COARSE FRAGMENTATION'
    c1 = 1; c2 = 1
    for line in psirr.readlines():
        if scfbas in line:
            check1 += c1*1
            c1 = -1; continue
        elif 'MAIN FRAGMENT' in line:
            check2 += c2*1
            c2 = -1; continue
        elif 'EUMP2' in line:
            check3 = 1
        elif 'SCF Done:  E(RHF)' in line:
            count2 += 1
            if chk == 0:
                if count2 == 1:
                    count2 = 0
                    check4 = 1
            elif chk == 1:
                if not param.coarse_level == 2:
                    check4 = 1
                else:
                    if count2 == 2:
                        check4 = 1
            elif chk == 2:
                if count2 == 1:
                    check4 = 1

        if check3 == 1:
            s = line.split()
            energymp2p = float(s[-1].replace('D','e'))
            check3 = 0
        elif check4 == 1:
            s = line.split()
            energyscfp = float(s[4])
            check4 = 0

        if check1 == 1:
            if 'Forces (Hartrees/Bohr)' in line:
                count3 = 1
                continue
            elif not count3:
                continue
            if (count3 and 'Number     Number              X              Y              Z' in line) or\
               (count3 and '------------------' in line):
                continue
            if count1+1 > len(f1.coord):
                check1 = 0
                count1 = 0
                count3 = 0
                continue
            else:
                s = line.split()
                s1 = [-float(xx) for xx in s]
                Grad.addc(f1.coord[count1].n,s1[2:])
                count1 += 1

        elif check2 == 1:
            if 'Forces (Hartrees/Bohr)' in line:
                count3 = 1
                continue
            elif not count3:
                continue
            if (count3 and 'Number     Number              X              Y              Z' in line)or\
               (count3 and '------------------' in line):
                continue
            if count1+1 > len(f1.coord):
                check2 = 0
                count1 = 0
                continue
            else:
                s = line.split()
                s1 = [-float(xx) for xx in s]
                Grad.addf(f1.coord[count1].n,s1[2:])
                count1 += 1

    psirr.close()
    if chk == 0:
        return(energyscfp,energymp2p,Grad.c,Grad.f)
    elif chk == 1 or chk == 2:
        return(energyscfp,Grad.c)

def s_gau(file1):
    check1 = 0
    check2 = 0
    check3 = 0
    check4 = 0
    check5 = 0
    count1 = 0
    count2 = 0
    Grad = Gradobj()
    Xcount = 0
    psirr = open(file1,'r')
    for line in psirr.readlines():
        if 'Forces (Hartrees/Bohr)' in line:
            check1 = 1
            continue
        elif 'SCF Done:  E(RHF)' in line:
            s = line.split()
            energyscfp = float(s[4])

        if check1 == 1:
            count1 += 1
            if count1 <=2 :
                continue
            if '----------' in line:
                check1 = 0
                continue
            s = line.split()
            s1 = [-float(xx) for xx in s]
            Grad.adds(0,s1[2:])                
    return(Grad.s,energyscfp)

def smp2_gau(file1):
    check1 = 0
    check2 = 0
    check3 = 0
    check4 = 0
    check5 = 0
    count1 = 0
    count2 = 0
    Grad = []
    Xcount = 0
    cout = 1
    psirr = open(file1,'r')
    for line in psirr.readlines():
        if 'Forces (Hartrees/Bohr)' in line:
            check1 = 1
            continue
        elif 'EUMP2' in line:
            s = line.split()
            energyscfm = float(s[-1].replace('D','e'))
        elif 'SCF Done:  E(RHF)' in line:
            s = line.split()
            energyscfp = float(s[4])

        if check1 == 1:
            count1 += 1
            if count1 <=2:
                continue
            if '----------' in line:
                check1 = 0
                continue
            s = line.split()
            s1 = [-float(xx) for xx in s]
            Grad.append(Gmethod(cout,s1[2:]))
            cout += 1
                
    return(Grad,energyscfm)

#-----------------------------------------------------
#-----------------------------------------------------

class FRAGMENT:
    def __init__(self,C):
        self.bon = Fmethod(C)
        self.nbon = Fmethod(C)
        if C == 'E':
            self.sup1 = 0.0
            if param.coarse:
                if param.coarse_level == 2:
                    self.sup2 = 0.0
        elif C == 'G':
            self.sup1 = []
            if param.coarse:
                if param.coarse_level == 2:
                    self.sup2 = []
            
    def addfB(self,E):
        self.bon.addf(E)

    def addcB(self,E):
        self.bon.addc(E)

    def addccB(self,E):
        self.bon.addcc(E)
        
    def addccaB(self,E):
        self.bon.addcca(E)
        
    def addfNB(self,E):
        self.nbon.addf(E)

    def addcNB(self,E):
        self.nbon.addc(E)

    def addccNB(self,E):
        self.nbon.addcc(E)
        
    def addccaNB(self,E):
        self.nbon.addcca(E)
        
    def correct1(self):
        self.sup1 = self.bon.cc + self.nbon.cc

    def correct2(self):
        self.cca = self.sup2 -(self.bon.cca + self.nbon.cca)

    def printfE(self,out):
        out.write('\n')
        out.write(' Fragment Energies \n')
        out.write(' Bonded contribution      : {:>13.8f} H\n'.format(self.bon.f))
        if param.con == 'fnb':
            out.write(' Non-bonded contribution  : {:>13.8f} H\n'.format(self.nbon.f))
        out.write(' Correction increment (1) : {:>13.8f} H\n'.format(self.sup1-
                                                        self.bon.c-self.nbon.c))
        if param.coarse:
            if param.coarse_level == 2:
                out.write(' Correction increment (2) : {:>13.8f} H\n'.format(self.cca))

    def Energy(self):
        E = self.sup1 + self.bon.f - self.bon.c

        E = self.bon.f + self.nbon.f + (self.sup1-self.bon.c-self.nbon.c)
        if param.coarse:
            if param.coarse_level == 2:
                E += self.cca
        return E
    
    def fragmentationcc(self):
        list1 = bonfrag(self.bon.cc,1); list2 = nbonfrag(self.nbon.cc,0)
        for i in range(param.natom):
            xg = list1[i].gx + list2[i].gx; yg = list1[i].gy + list2[i].gy
            zg = list1[i].gz + list2[i].gz
            self.sup1.append(Gmethod(i+1,[xg,yg,zg]))
    
    def fragmentationcca(self):
        L1 = bonfrag(self.bon.bon.cca,1); L2 = nbonfrag(self.nbon.cca,0)
        emptylist = []
        for i in range(param.natom):
            xg = self.sup2[i].gx - (L1[i].gx + L2[i].gx) 
            yg = self.sup2[i].gy - (L1[i].gy + L2[i].gy)
            zg = self.sup2[i].gz - (L1[i].gz + L2[i].gz)
            emptylist.append(Gmethod(i+1,[xg,yg,zg]))
        return emptylist

class Fmethod:
    def __init__(self,C):
        self.add = False; self.append = False
        if C == 'E':
            self.add = True
        elif C == 'G':
            self.append = True
        
        if self.add:
            self.f = 0.0
            self.c = 0.0
            if param.coarse:
                self.cc = 0.0
                if param.coarse_level == 2:
                    self.cca = 0.0
        elif self.append:
            self.f = []
            self.c = []
            if param.coarse:
                self.cc = []
                if param.coarse_level == 2:
                    self.cca = []

    def addf(self,E):
        if self.add:
            self.f += E
        elif self.append:
            self.f.append(E)

    def addc(self,E):
        if self.add:
            self.c += E
        elif self.append:
            self.c.append(E)
        
    def addcc(self,E):
        if self.add:
            self.cc += E
        elif self.append:
            self.cc.append(E)

    def addcca(self,E):
        if self.add:
            self.cca += E
        elif self.append:
            self.cca.append(E)

class Gradobj:
    
    def __init__(self):
        self.c = []
        self.f = []
        self.s = []

    def addc(self,n,g):
        self.c.append(Gmethod(n,g))
    
    def addf(self,n,g):
        self.f.append(Gmethod(n,g))

    def adds(self,n,g):
        self.s.append(Gmethod(n,g))

class Gmethod:
    
    def __init__(self,n,g):
        self.n = n
        self.gx = g[0]
        self.gy = g[1]
        self.gz = g[2]
        


# ---> Notes



'''
NOTE1: 
1. Neglecting non-bonded contributions in the correction term can lead
   to an optimization without any interaction between the framents. 
   This is risky since the fragments can fall on each other resulting
   in high energy.
2. This was tested for many cases of folded poly peptide chain. See the
   supporting info of 
   Journal of Theoretical and Computational Chemistry 17 (05), 1850037
'''
