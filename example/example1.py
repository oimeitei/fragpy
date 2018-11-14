#!/usr/bin/env python
from fragpy.min import optf
from fragpy.opt_param import param

'''  RUN THIS PYHON SCRIPT 

1. Supply the geometry in the cartesian coordinates as
   a list in the format shown below.
   --------------
   READ the second point in exampl1.inp
   The format of the example given there would be
   coord = [ 
    ['C1',  0.0,    0.0,    0.0,   1],
    ['C2',  0.0,    0.0,    0.0,   1],
    ['C3',  0.0,    0.0,    0.0,   2],
    ['C4',  0.0,    0.0,    0.0,   2],
    ['C5',  0.0,    0.0,    0.0,   3],
    ['C6',  0.0,    0.0,    0.0,   3]]
    -------------

2. Set the keywords from KEYLIST to overwrite the 
   default parameters 
        param.KEY = key

3. call optf(geometry list) as shown below

'''

coord = [
['N1', -1.08192199,  7.43133438, -3.64184138],
['H1', -0.79076004,  7.19359607, -2.68385020],
['H1', -0.49113179,  6.93510001, -4.32235455],
['C1', -2.49213835,  7.20527042, -3.86446890],
['H1', -3.08131841,  7.84711944, -3.17857026],
['H1', -2.75433466,  7.53009208, -4.89260723],
['O1', -4.12999967,  5.43099350, -3.71446813],
['C1', -2.94499239,  5.76007962, -3.68429206],
['N1', -1.94917148,  4.81923883, -3.51666770],
['H1', -1.00481229,  5.12920805, -3.22129286],
['C1', -2.29382867,  3.48607934, -3.10329606],
['H1', -3.29923581,  3.48775881, -2.62662227],
['H1', -2.36419144,  2.78421120, -3.97189471],
['O2', -0.24550297,  3.69038154, -1.86045264],
['C2', -1.22708406,  3.00017838, -2.14324166],
['N2', -1.39516610,  1.73976680, -1.62557414],
['H2', -2.29182440,  1.23023370, -1.71858522],
['C2', -0.51292429,  1.24903257, -0.60302627],
['H2',  0.44856117,  0.87383198, -1.03236614],
['H2', -0.23488799,  2.06631712,  0.10011707],
['O2', -2.35472921, -0.23954767, -0.19233017],
['C2', -1.22398364,  0.12695505,  0.12552513],
['N2', -0.54481297, -0.45131146,  1.17139420],
['H2',  0.46612149, -0.28359419,  1.31595027],
['C2', -1.07654954, -1.62047446,  1.81665562],
['H2', -1.81101489, -1.35393145,  2.61617210],
['H2', -1.63528965, -2.24430432,  1.08341291],
['O3',  1.24485827, -1.99446223,  2.31636679],
['C3',  0.07616160, -2.39324742,  2.41857561],
['N3', -0.23946707, -3.53968896,  3.09519996],
['H3', -1.20161960, -3.90278700,  3.07849861],
['C3',  0.79476081, -4.35934650,  3.69797601],
['H3',  1.52347849, -3.69611439,  4.20739718],
['H3',  0.32969871, -5.01695306,  4.45201078],
['O3',  1.32430557, -6.40775533,  2.51197580],
['C3',  1.50646153, -5.19728098,  2.63675865],
['N3',  2.38561170, -4.51965004,  1.83371575],
['H3',  2.37132330, -3.48157744,  1.85154108],
['C3',  2.91895020, -5.19994650,  0.67737908],
['H3',  2.31956534, -4.97337818, -0.23520590],
['H3',  3.95333916, -4.84124263,  0.48087434],
['C3',  2.96591068, -6.70953428,  0.91523571],
['H3',  3.58250126, -7.02900092,  1.82807634],
['O3',  2.56279999, -7.49900720,  0.10693198]]


param.ITMAX = 50
optf(coord)
