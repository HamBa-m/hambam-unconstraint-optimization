from compare import comp_unidim,comp_multivar,comp_sys_solv
from functions import g,gpp,gp,h,H

"""
the numbers of the 10 operations for unidimensional optimization:
1: usfs 
2: usas
3: bf
4: bs
5: halv
6: Fibo
7: golden
8: Newton-Rapson
9: Quasi-Newton
10: Secant
"""
comp_unidim(g,-3,5,[0,1,2,3])
comp_unidim(g,-2,4,[4,5,6])
comp_unidim(g,-7,1,[7,8,9],gp,gpp)

comp_multivar(h, [3.,5.],H,10)

A,B = [[5,2],[2,3]],[1,1]
comp_sys_solv(A,B)