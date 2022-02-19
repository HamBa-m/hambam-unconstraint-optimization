from unidim import elimination as e
from unidim import interpolation as i
from multivar import gauss as gs
from multivar import deco as dc
from multivar import gradient as gd
from time import perf_counter
import numpy as np
from matplotlib import pyplot as plt

func = [e.usfs,e.usas,e.bf,e.bs,e.halv,e.fibo,e.golden,i.Newton_Rapson,i.Quasi_Newton,i.Secant_Method,]
func2 = [gd.gradient_descent,gd.conjugate_gradient,gd.newton_descent]
func3 = [gs.SolveMatSys,dc.solveLU,dc.solveCholesky]

def comp_unidim(f,a,b,X,fp=None,fpp=None):
    """
the numbers associated to the 10 operations for unidimensional optimization:
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
    T = list()
    for u in X:
        if u < 7 or u == 8: 
            t = perf_counter()
            func[u](f,a,b)
            T.append(perf_counter() - t)
        elif u == 7 :
            t = perf_counter()
            func[u](f,fp,fpp,a,b)
            T.append(perf_counter() - t)
        else :
            t = perf_counter()
            func[u](f,fp,a,b)
            T.append(perf_counter() - t)
    i = 0
    for e in T:
        plt.barh(i,e, 0.3, label = str(func[X[i]]))
        i += 1
    plt.xlabel("time in (s)")
    plt.legend()
    plt.grid()
    plt.show()
    return

def comp_multivar(f, x, F, axlim):
    T = list()
    for u in range(3):
        t = perf_counter()
        func2[u](f,x,F=F,axlim=axlim)
        T.append(perf_counter() - t)
    i = 0
    for e in T:
        plt.barh(i,e, 0.3, label = str(func2[i]))
        i += 1
    plt.xlabel("time in (s)")
    plt.legend()
    plt.grid()
    plt.show()
    return


def comp_sys_solv(A,B):
    T = list()
    for u in range(3):
        t = perf_counter()
        func3[u](A,B)
        T.append(perf_counter() - t)
    i = 0
    for e in T:
        plt.barh(i,e, 0.3, label = str(func3[i]))
        i += 1
    plt.xlabel("time in (s)")
    plt.legend()
    plt.grid()
    plt.show()
    return 