from time import perf_counter
from math import sqrt, ceil, log
from matplotlib import pyplot as plt
import numpy as np

#unrestricted search with fixed step
def usfs(f, a, b, p = 0.01, plot=False):
    """
    f: monovariable function
    a : start of search interval
    b : end of search interval
    p : step size
    plot : a parameter to choose either to plot or not the results
    """
    if plot :
        X = np.linspace(a,b,1000)
        plt.plot(X, f(X))
    t = perf_counter()
    iter = 0
    A = []
    while f(a) > f(a + p) and a <= b :
        iter += 1
        a += p
        if plot and len(A) < 10 : A.append(a)
    t = perf_counter() - t
    print("time execution: ",t)
    if plot :
        for e in A :
            plt.plot(e,f(e),"o")
    if a > b : 
        plt.plot(a-p, f(a-p),"s",label="optimum")
        plt.show()
        return a - p
    if plot : 
        plt.plot(a,f(a),"s",label="optimum")
        plt.title("optimization with fixed step search")
        plt.legend()
        plt.show()
    return a

#unrestricted search with accelerated step
def usas(f, a, b, p = 0.01, plot=False):
    """
    f: monovariable function
    a : start of search interval
    b : end of search interval
    p : step size
    plot : a parameter to choose either to plot or not the results
    """
    t = perf_counter()
    it , step = a, 0
    A = []
    while  step != p :
        if it > b : it  -= step / 2
        step , i1, i2 = p, f(it), f(it + p)
        while i1 > i2 and it < b :
            if plot and len(A) < 10 : A.append(it)
            it += step
            i1 = i2
            step *= 2
            i2 = f(it + step)
    t = perf_counter() - t
    print("time execution: ",t)
    if plot :
        plt.plot(it,i1,"s",label="optimum")
        X = np.linspace(a,b,1000)
        plt.plot(X,f(X))
        for e in A :
            plt.plot(e,f(e),"o")
        plt.legend()
        plt.title("optimization with accelerated step search")
        plt.show()
    return it

#exhaustive search
def bf(f, a, b, n = 1000, plot=False):
    """
    f: monovariable function
    a : start of search interval
    b : end of search interval
    n : number of divisions of the interval
    plot : a parameter to choose either to plot or not the results
    """
    t = perf_counter()
    p = (b - a) / (n - 1)
    it = a
    L = []
    A = []
    while a <= b :
        L.append(f(a))
        a += p
        if plot : A.append(a)
    t = perf_counter() - t
    print("time execution: ",t)
    if plot : 
        plt.plot(it + L.index(min(L))*p,f(it + L.index(min(L))*p),"s",label="optimum")
        X = np.linspace(a,b,1000)
        plt.plot(X,f(X))
        for e in A :
            plt.plot(e,f(e),"o")
        plt.legend()
        plt.title("optimization with Brute Force method")
        plt.show()
    return it + L.index(min(L))*p

#binary search
def bs(f, a, b, e = 0.01, delta = 0.001, plot=False):
    """
    f: monovariable function
    a : start of search interval
    b : end of search interval
    e : search precision
    delta : a factor relied to the binary search
    plot : a parameter to choose either to plot or not the results
    """
    t = perf_counter()
    A = []
    mid = (a + b) / 2
    n = ceil(log((b - a - delta) / (2 * e * (b - a) - delta)) * (2 / log(2)))
    if n%2 : n+=1
    if plot : 
        X = np.linspace(a,b,1000)
        plt.plot(X,f(X))
    for i in range(n) :
        y1, y2 = f(mid - (delta / 2)), f(mid + (delta / 2))
        if y2 > y1 : b = mid + (delta / 2)
        else : a = mid - (delta / 2)
        mid = (a + b) / 2
        if plot : A.append(mid)
    t = perf_counter() - t
    print("time execution: ",t)
    if plot : 
        plt.plot(mid,f(mid),"s",label="optimum")
        for u in A :
            plt.plot(u,f(u),"o")
        plt.legend()
        plt.title("optimization with Binary Search method")
        plt.show()
    return mid

#interval halving method
def halv(f, a, b, e = 0.01, plot=False):
    """
    f: monovariable function
    a : start of search interval
    b : end of search interval
    e : search precision
    plot : a parameter to choose either to plot or not the results
    """
    t = perf_counter()
    A = [[],[],[]]
    n = ceil((2 * log(2 * e))/ log(1 / 2)) + 1
    if (n%2)==0 : n += 1
    if plot :
        X = np.linspace(a,b,1000)
        plt.plot(X,f(X))
    for i in range(n):
        x0, x1, x2 = (a+b)/2, (a+b)/3, (2*(a+b))/3
        y0, y1, y2 = f(x0), f(x1), f(x2)
        if plot :
            A[0].append(x0)
            A[1].append(x1)
            A[2].append(x2)
        if y1 < y0 < y2 : b = x0
        elif y1 > y0 > y2 :	a = x0 
        else :
            a = x1
            b = x2
    t = perf_counter() - t
    print("time execution: ",t)
    if plot :
        plt.plot((a+b)/2, f((a+b)/2),"s",color="yellow",label="optimum")
        for i in range(n):
            plt.plot(A[0][i],f(A[0][i]),"o",color='red')
            plt.plot(A[1][i],f(A[1][i]),"o",color='blue')
            plt.plot(A[2][i],f(A[2][i]),"o",color='green')
        plt.legend()
        plt.title("optimization with Interval Halving method")
        plt.show()
    return (a + b) / 2

#fibonacci method
def fibo(f, a, b, n = 20, plot=False):
    """
    f: monovariable function
    a : start of search interval
    b : end of search interval
    n : number of iterations (deepth)
    plot : a parameter to choose either to plot or not the results
    """
    A = [[],[]]
    if plot :
        X = np.linspace(a,b,1000)
        plt.plot(X,f(X))
    t = perf_counter()
    f0, f1, f2 = 1, 1, 2
    for i in range(n-2) :
        f0 = f1
        f1 = f2
        f2 = f1 + f0
    c = a + (f0 / f2) * (b - a)
    d = a + (f1 / f2) * (b - a)
    y1 = f(c)
    y2 = f(d)
    for i in range(n):
        if plot and len(A[0]) < 10 : 
            A[0].append(c)
            A[1].append(d)
        if y1 < y2 :
            b = d
            d = c
            c = a + (b - d)
            y2 = y1
            y1 = f(c)
        else :
            a = c
            c = d
            d = a + (b - c)
            y1 = y2
            y2 = f(d)
    t = perf_counter() - t
    print("time execution: ",t)
    if plot :
        plt.plot(c,f(c),"s",label="optimum")
        for i in range(len(A[0])):
            plt.plot(A[0][i],f(A[0][i]),"o",color='red')
            plt.plot(A[1][i],f(A[1][i]),"o",color='blue')
        plt.legend()
        plt.title("optimization with Fibonacci method")
        plt.show()
    return c

#golden ratio method
def golden(f, a, b, e = 0.01, plot=False):
    """
    f: monovariable function
    a : start of search interval
    b : end of search interval
    e : search precision
    plot : a parameter to choose either to plot or not the results
    """
    A = [[],[]]
    if plot and len(A[0]) < 10 :
        X = np.linspace(a,b,1000)
        plt.plot(X,f(X))
    t = perf_counter()
    phi = (sqrt(5) - 1) / 2
    Ln = 1
    f0, f1, f2 = 1, 1, 2
    for i in range(100) :
        f0 = f1
        f1 = f2
        f2 = f1 + f0
    c = a + (f0 / f2) * (b - a)
    d = a + (f1 / f2) * (b - a)
    y1 = f(c)
    y2 = f(d)
    while Ln > e:
        if plot : 
            A[0].append(c)
            A[1].append(d)
        if y1 < y2 :
            b = d
            d = c
            c = a + (b - d)
            y2 = y1
            y1 = f(c)
        else :
            a = c
            c = d
            d = a + (b - c)
            y1 = y2
            y2 = f(d)
        Ln *= phi
    t = perf_counter() - t
    print("time execution: ",t)
    if plot :
        plt.plot(c,f(c),"s",label="optimum")
        for i in range(len(A[0])):
            plt.plot(A[0][i],f(A[0][i]),"o",color='red')
            plt.plot(A[1][i],f(A[1][i]),"o",color='blue')
        plt.title("optimization with Golden section")
        plt.legend()
        plt.show()
    return c
