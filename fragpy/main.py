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
import os,sys,numpy,pickle,pxssh,getpass,glob,re
from fragpy.opt_param import param
from fragpy.write import splitwritef
from fragpy.grepfiles import fraggradpsi
from fragpy.out_print import out_error
from subprocess import Popen,PIPE
from pexpect import *
from multiprocessing import Pool,Process,Pipe,Queue,Manager,Value,Array,Lock 



def mainfunc(Q):
    # Write inputs
    inpfil,frags = splitwritef(Q) 
 
    # Run
    if param.con == 's':
        runpsi(inpfil)
        apple = 0
    elif param.con == 'fnb':
        runpsi(inpfil)
        apple = 0
    elif param.con == 'fb':
        runpsi(inpfil)
        apple = 0

    # Read & get gradients, E
    glist,etot = fraggradpsi(frags)
    return(glist,etot)

def run1(r1cn,r1lck,r1node,inputfil):
    pathp=os.getcwd()
    r1 = pxssh.pxssh(timeout=None)
    r1.force_password = True
    r1.login(r1node,param.USER)
    r1.logfile=open('logo','w')
    r1.sendline("cd "+str(pathp))
    r1.prompt()

    for bl in param.bashlines:
         r1.sendline(bl)
         r1.prompt()

    if param.program.upper() == 'PSI4': 
         scvaria1 = param.scvaria+'/'+str(r1node)
         if not os.path.isdir(scvaria1):
              os.mkdir(scvaria1)
         r1.sendline("export PSI_SCRATCH="+scvaria1)
         r1.prompt()

    r1conv1 = False
    while not r1conv1:
        r1lck.acquire()
        r1cn.value += 1
        valr1 = r1cn.value - 1
        r1lck.release()
        if param.program.upper() == 'PSI4':
            r1.sendline(param.PSI4+" -i "+inputfil[valr1]+" -o "+
                         inputfil[valr1].split('.')[0]+".out -n "+str(param.procs))
        
        elif param.program.upper() == 'MOLPRO':
            r1.sendline(param.MOLPRO+' -t '+str(param.procs)+' '+inputfil[valr1])
        elif param.program.upper() == 'GAUSSIAN':
            r1.sendline(param.GAUSSIAN+' '+inputfil[valr1])
        else:   
             out_error('ERROR: No program to run Fragpy - set param.program in opt_param.py\n') 

        r1.prompt(timeout=None)
        r1lck.acquire()
        valr1 = r1cn.value
        r1lck.release()
        if valr1 == len(inputfil):
            r1conv1 = True

def run2(r2cn,r2lck,inputfil):
    r2conv2 = False
    while not r2conv2:
        r2lck.acquire()
        r2cn.value += 1
        valr2 = r2cn.value - 1
        r2lck.release()
        if param.program.upper() == 'PSI4':
            rrr = Popen([param.PSI4+' -i '+inputfil[valr2]+' -o '+inputfil[valr2].
                          split('.')[0]+'.out -n '+\
                          str(param.procs)],shell=True)
        
        elif param.program.upper() == 'MOLPRO':
            rrr = Popen([param.MOLPRO+' -t '+str(param.procs)+' '+inputfil[valr2]],shell=True)
        elif param.program.upper() == 'GAUSSIAN':
            rrr = Popen([param.GAUSSIAN+' '+inputfil[valr2]],shell=True)
        else:
            out_error('ERROR: No program to run Fragpy - set param.program in opt_param.py\n') 
        
        rrr.wait()
        r2lck.acquire()
        valr2 = r2cn.value
        r2lck.release()
        if valr2 == len(inputfil):
            r2conv2 = True

#----Hybrid parallelisation
# omp on each compute node spread over 'nodes'
# the only communication reqd. is rnmcounter & rnmlock
# run2 on the first allocated node
# run1 on the rest

def runpsi(inputfil): 

     if param.parallel > 1:
         nodefil = open('nodes.txt','r')
         nodes = []
         for i in nodefil.readlines():
              s = i.split()
              for j in s:
                   if not j in nodes:
                        nodes.append(j)

     rnm = Manager()
     rnmcounter = rnm.Value('i',0) # index of inputfil
     rnmlock = rnm.Lock()
     pool1 = Pool()
     if param.parallel > 1:
         pool2 = []
         len_nodes = len(nodes)-1
         for rpsi in range(len_nodes):
              pool2.append(Pool())
         rs2 = []
         for rpsi in range(len_nodes):
              rs2.append(pool2[rpsi].apply_async(run1,[rnmcounter,rnmlock,nodes[rpsi+1],inputfil]))
     
     rs1=pool1.apply_async(run2,[rnmcounter,rnmlock,inputfil])


     if param.parallel > 1:
          for rpsi in rs2:
               rpsi.get()
     rs1.get()

     if param.parallel > 1:
          for rpsi in pool2:
               rpsi.close()
     pool1.close()

