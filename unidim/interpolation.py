from time import perf_counter
import numpy as np
from matplotlib import pyplot as plt

#returns the minimum of a unidimensional unimodal function using Newton Rapson method
def Newton_Rapson(f,fp, fpp, a, e = 0.1, plot=False):
    t = perf_counter()
    A = []
    k = a
    while abs(fp(k)) > e :
        k = k - fp(k)/fpp(k)
        if plot :
            A.append(k)
    t = perf_counter() - t
    print("time execution: ", t)
    if plot : 
        X = np.linspace(a,2*k-a,1000)
        plt.plot(X,f(X))
        plt.plot(k,f(k),"s",label="optimum")
        for u in A:
            plt.plot(u,f(u),"o")
        plt.legend()
        plt.show()
    return k

#returns the minimum of a unidimensional unimodal function using Quasi Newton method
def Quasi_Newton(f, a, s = 0.01, e = 0.1, plot=False):
    t = perf_counter()
    A = []
    k = a
    dev = abs((f(k + s) - f(k - s))/(2 * s))
    while  dev > e:
        if plot :
            A.append(k)
        x = f(k + s)
        y = f(k - s)
        z = f(k)
        k = k - ((s * (x - y))/(2 * (x - 2 * z + y)))
        dev = abs((f(k + s) - f(k - s))/(2 * s))
    t = perf_counter() - t
    print("time execution: ", t)
    if plot : 
        X = np.linspace(a,2*k-a,1000)
        plt.plot(X,f(X))
        for u in A:
            plt.plot(u,f(u),"o")
        plt.plot(k,f(k),"s",label="optimum")
        plt.legend()
        plt.show()
    return k

#returns the minimum of a unidimensional function usingthe Secant method
def Secant_Method(f, fp, lamda, s = 0.01, e = 0.01, plot=False):
    A = []
    t = perf_counter()
    a = lamda
    b = a + s
    while fp(s) < 0:
        a = s
        s *= 2
    b = s
    k = a - ((fp(a) * (b - a)) / (fp(b) - fp(a)))
    while abs(fp(k)) > e:
        if plot :
            A.append(k)
        if fp(k) >= 0 :
            b = k
        else : a = k
        k = a - ((fp(a) * (b - a)) / (fp(b) - fp(a)))
    t = perf_counter() - t
    print("time execution: ", t)
    if plot : 
        X = np.linspace(lamda,2*k-lamda,1000)
        plt.plot(X,f(X))
        for u in A:
            plt.plot(u,f(u),"o")
        plt.plot(k,f(k),"s",label="optimum")
        plt.legend()
        plt.show()
    return k