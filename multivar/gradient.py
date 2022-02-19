import numpy as np 
import scipy.optimize as op 
import scipy.linalg as la
import numdifftools as nd 
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
from time import perf_counter

#returns the minimum of a multivariable function using gradient descent
def gradient_descent(f, xk, delta = 0.01, plot=False, F = None, axlim = 10):
    if plot : ax = plt.axes(projection='3d')
    A = []
    t = perf_counter()
    dk = nd.Gradient(f)(xk)
    while la.norm(dk) > delta :
        if plot and len(A) < 10 : A.append(xk)
        xt = xk
        phi = lambda s : f(xk - s * dk)
        alpha = op.newton(phi, 1)
        xk -= alpha * dk
        if plot and len(A) < 10 : A.append(xk)
        dk = nd.Gradient(f)(xk)
        if la.norm(xk - xt) < delta : break
    t = perf_counter() - t
    print("execution time: ",t)
    if plot :
        for u in A:
            ax.scatter(u[0], u[1], f(u), c = 'b', s = 50)
        ax.scatter(xk[0], xk[1], f(xk), c = 'r', s = 50,label="optimum")
        x = np.arange(-axlim, axlim, axlim/100)
        y = np.arange(-axlim, axlim, axlim/100)
        X, Y = np.meshgrid(x, y)
        Z = F(X,Y)
        ax.set_xlabel('x', labelpad=20)
        ax.set_ylabel('y', labelpad=20)
        ax.set_zlabel('z', labelpad=20)
        surf = ax.plot_surface(X, Y, Z, cmap = plt.cm.cividis)
        plt.legend()
        plt.title("optimizition with Gradient Descent")
        plt.show()
    return xk

#returns the minimum of a multivariable function using conjugate gradient
def conjugate_gradient(f, x, plot=False, F = None,axlim = 10):
    if plot : ax = plt.axes(projection='3d')
    A = []
    t = perf_counter()
    d = -nd.Gradient(f)(x)
    q = nd.Hessian(f)(x)
    n = len(x)
    for k in range(1, n):
        if plot and len(A) < int(n/3) : A.append(x)
        alpha = (d.T@d)/(d.T@q@d)
        x += alpha * d
        beta = (nd.Gradient(f)(x).T@q@d)/(d.T@q@d)
        d = beta * d - nd.Gradient(f)(x)
    t = perf_counter() - t
    print("execution time: ",t)
    if plot :
        for u in A:
            ax.scatter(u[0], u[1], f(u), c = 'b', s = 50)
        ax.scatter(x[0], x[1], f(x), c = 'r', s = 50,label="optimum")
        x = np.arange(-axlim, axlim, axlim/100)
        y = np.arange(-axlim, axlim, axlim/100)
        X, Y = np.meshgrid(x, y)
        Z = F(X,Y)
        ax.set_xlabel('x', labelpad=20)
        ax.set_ylabel('y', labelpad=20)
        ax.set_zlabel('z', labelpad=20)
        surf = ax.plot_surface(X, Y, Z, cmap = plt.cm.cividis)
        plt.legend()
        plt.title("optimizition with Conjugate Gradient")
        plt.show()
    return x

#verifies if a function is defined positive at a point x
def is_pos_def(f, x):
    m = nd.Hessian(f)(x)
    return np.all(np.linalg.eigvals(m) > 0)

#returns the minimum of a multivariable function using Newton descent
def newton_descent(f, x, delta = 0.01, plot=False, F = None, axlim = 10):
    d = -la.inv(nd.Hessian(f)(x))@nd.Gradient(f)(x)
    m = nd.Hessian(f)(x)
    if la.det(m) == 0:
        return ValueError
    else :
        if plot : ax = plt.axes(projection='3d')
        A = []
        t = perf_counter()
        while la.norm(d) > delta:
            if plot and len(A) < 10 : A.append(x)
            phi = lambda s : f(x - s * d)
            alpha = op.newton(phi, 1)
            x -= alpha * d
            dt = nd.Hessian(f)(x)
            if is_pos_def(f, x) : d = -la.inv(nd.Hessian(f)(x))@nd.Gradient(f)(x)
            else :
                u = la.eigvals(f)(x)
                epsilon = min(u)
                n = len(x)
                d = -la.inv((epsilon * np.identity(n) + nd.Hessian(f)(x))) @ nd.Gradient(f)(x)
        t = perf_counter() - t
        print("execution time: ",t)
        if plot :
            for u in A:
                ax.scatter(u[0], u[1], f(u), c = 'b', s = 50)
            ax.scatter(x[0], x[1], f(x), c = 'r', s = 50,label="optimum")
            x = np.arange(-axlim, axlim, axlim/100)
            y = np.arange(-axlim, axlim, axlim/100)
            X, Y = np.meshgrid(x, y)
            Z = F(X,Y)
            ax.set_xlabel('x', labelpad=20)
            ax.set_ylabel('y', labelpad=20)
            ax.set_zlabel('z', labelpad=20)
            surf = ax.plot_surface(X, Y, Z, cmap = plt.cm.cividis)
            plt.legend()
            plt.title("optimizition with Newton descent")
            plt.show()
    return x

