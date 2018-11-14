Fragpy
======
Fragpy is a tool to optimize geometries of large molecules using *Incremental Molecular Fragmentation* method of OR Meitei & A Hesselmann, *J. Theor. Comput. Chem.*, 17(05), 185007. Fragpy works with other quantum chemistry programs such as PSI4, Molpro, Gaussian. Written in Python, it can perform fragment calculations on several compute nodes with minimal interconnect between the nodes. Fragpy works with MP2 method (and HF method for the correction term) for the molecular fragmentation but it should be easy to adapt Fragpy for other methods as well.

## Installation
* Requirements
    - Python 2.7
    - Numpy 1.15.1
    - PSI4 / Molpro / Gaussian program.
* In fragpy/opt_param.py set the path of the executable (as a string) of PSI4/ Molpro/ Gaussian. 

         param.PSI4 = "/path/to/psi4/bin/psi4"
* In parallel execution of Fragpy over more than one nodes, bash line commands can be passed (which might become necessary) by including the command lines (as a string) as a list in fragpy/opt_param.py.
     	 
         param.bashlines = ['*some bash command*','*another bash command*']
* Just inclide the top level directory in enviromental variable "PYTHONPATH".
  
         export PYTHONPATH=/path/to/fragpy:PYTHONPATH	
## Usage
Fragpy is very flexible to use, see fragpy/example. Also read fragpy/KEYLIST for a list of keywords. 