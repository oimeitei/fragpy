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
Default parameters for the fragment geometry optimisation
'''
import getpass,os

class optf:
     def __init__(self):

          # ---> fragpy parameters         
          self.con = 'fnb'
          self.DEBUG = 1
          self.Prestart = False
          self.restart = False

          # ---> Method
          self.method = 'mp2'
          self.cmethod = 'scf'

          # ---> Basis
          self.basis = 'cc-pVDZ'
          self.basis_jk = 'cc-pVDZ-jkfit'
          self.basis_mp2 = 'cc-pVDZ-ri'
          self.cbasis = 'STO-3G'
          self.cbasis_jk = 'cc-pVDZ-jkfit'
          self.cbasis_mp2 = 'cc-pVDZ-ri'

          # ---> Convergence
          self.E_conv = '1e-8'   # SCF
          self.D_conv = '1e-8'   # SCF
          self.O_conv = '1e-7'   # SCF
          self.MAXG = 3.0e-4
          self.DELE = 1.0e-6
          self.MAXQ = 1.8e-3

          # ---> Cutoffs
          self.bonthres = 4.0e0     
          self.cbonthres = 4.0e0         
          self.thresnb  = 4.0e10    
          self.cthresnb = 4.0e10

          # ---> Opt parameters
          self.old_trust = False    
          self.OPTstep='RF'    
          self.modelHtype = 'R' 
          self.LBFGS = True
          self.ITMAX=300
          self.trust = 0.3
          self.stepmax = 0.3
          self.stepmin = 0.01
          self.updateT = False
          self.updateFT = False
          self.updateH = False

          # ---> Fragment parameters
          self.FIXGEOM = True
          self.l1o1 = True          
          self.nbscale = 1e0 
          self.coarse = False
          self.coarse_level = 2 
          self.dumfrag = False
          self.dumfrag_level = 1 # No. of adjoining groups
          
          # ---> System related
          self.program = 'PSI4'
          self.parallel = 1
          self.procs = 1
          self.USER = getpass.getuser()

          self.memory = 1000
          self.MOLPRO = ''
          self.PSI4 = ''
          self.GAUSSIAN = ''
          self.gau_command = []
          
          
param = optf()


param.scvaria = str(os.environ["SCVARIA"])
