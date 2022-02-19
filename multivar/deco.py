from time import perf_counter
from copy import deepcopy
from math import sqrt
from matrix import matrix as m
from multivar import gauss as gs

#returns the LU decomposition of a matrix
def lu(A):
    U = deepcopy(A)
    t = perf_counter()
    n = len(A)
    L = m.IdentityMat(n)
    for i in range(n):
        for j in range(i + 1, n):
            L[j][i] = U[j][i]/U[i][i]
            m.MatRowTransvection(U, j, i, -L[j][i])
    t = perf_counter() - t
    print("time execution for LU decomposition: ",t)
    return L, U

#returns true if the matrix is symetric
def is_symetric(A):
    for i in range(len(A)):
        for j in range(len(A)):
            if A[i][j] != A[j][i] : 
                return False
    return True

#returns true if the matrix is definite positive
def is_definite_postive(B):
    B = gs.GaussPivot(B)
    B = m.TransMat(B)
    B = gs.GaussPivot(B)
    for i in range(len(B)):
        if B[i][i] < 0 : return False
    return True

#returns the Cholesky's decomposition of a matrix
def choleski(A):
    if not(is_symetric(A)) or not(is_definite_postive(A)) :
        return ValueError, ValueError
    t = perf_counter()
    n = len(A)
    L = m.ZerosMat(n)
    for row in range(n):
        for col in range(row + 1):
            if (col == row):
                s = 0
                for i in range(col):
                    s += L[col][i]**2
                L[row][col] = sqrt(A[row][col] - s)
            else :
                s = 0
                for i in range(col):
                    s +=  L[row][i]*L[col][i]
                L[row][col] = (A[row][col] - s)/L[col][col]
    t = perf_counter() - t
    print("time execution for Cholesky decomposition: ",t)
    return L, m.TransMat(L)

#executes the forward substitution of a matrix
def ForSubst(L, B):
    n = len(B)
    Y = [0]*n
    Y[0] = B[0]/L[0][0]
    for i in range(1, n):
        s = 0
        for k in range(i):
            s += L[i][k] * Y[k]
        Y[i] = (B[i] - s)/L[i][i]
    return Y

#execute the backward substitution of a matrix
def BacSubst(U, Y):
    n = len(Y)
    X = [0]*n
    X[n - 1] = Y[n - 1]/U[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        s = 0
        for k in range(i + 1, n):
            s += U[i][k] * X[k]
        X[i] = (Y[i] - s)/U[i][i]
    return X

#solves a linear system using LU decomposition
def solveLU(A,B):
    t = perf_counter()
    L, U = lu(A)
    Y = ForSubst(L, B)
    Z = BacSubst(U, Y)
    t = perf_counter() - t
    print("time execution for solving the linear system using LU: ",t)
    return Z

#solves a linear system using Cholesky's decomposition
def solveCholesky(A,B):
    t = perf_counter()
    L, U = choleski(A)
    Y = ForSubst(L, B)
    Z = BacSubst(U, Y)
    t = perf_counter() - t
    print("time execution for solving the linear system using Cholesky: ",t)
    return Z