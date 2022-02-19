from copy import deepcopy
from time import perf_counter
from matrix import matrix as ma

#solving linear system with Gauss-Jordan
def SolveMatSys(N , C):
    M = deepcopy(N)
    B = deepcopy(C)
    t = perf_counter()
    row = 0
    col = 0
    pivot = 0
    n = len(M)
    m = len(M[0])
    while (row < n) and (col < m):
        if M[row][col] != 0 :
            temp = M[row][col]
            A = []
            for j in range(row + 1, n):
                A.append(M[j][col])
            for j in range(row + 1, n):
                ma.MatRowDilatation(M, j, temp)
                B[j] *= temp
            for j in range(row + 1, n):
                ma.MatRowTransvection(M, j, row, (-A[j - row - 1]))
                B[j] += B[row] * (-A[j - row - 1])
            row += 1
            col += 1
        else :
            ind = False
            for j in range(row + 1, n):
                if M[j][col] != 0:
                    ma.MatRowSwap(M, j, row)
                    ma.MatRowSwap(B, j, row)
                    ind = True
                    break
            if ind == False :
                col += 1
                row += 1
    row = len(B) - 1
    col = len(M[0]) - 1
    while row >= 0 and col >= 0 :
        if M[row][col] != 1 :
            B[row] /= M[row][col]
            M[row][col] = 1
        for i in range(row - 1, -1, -1):
            B[i] -= B[row] * M[i][col]
            M[i][col] = 0
        col -= 1
        row -= 1
    t = perf_counter() - t
    print("time execution: ",t)
    return B

#trigonalizing a matrice with Gauss-Jordan
def GaussPivot(N):
    M = deepcopy(N)
    t = perf_counter()
    row = 0
    col = 0
    pivot = 0
    n = len(M)
    m = len(M[0])
    while (row < n) and (col < m):
        if M[row][col] != 0 :
            temp = M[row][col]
            A = []
            for j in range(row + 1, n):
                A.append(M[j][col])
            for j in range(row + 1, n):
                ma.MatRowDilatation(M, j, temp)
            for j in range(row + 1, n):
                ma.MatRowTransvection(M, j, row, (-A[j - row - 1]))
            row += 1
            col += 1
        else :
            ind = False
            for j in range(row + 1, n):
                if M[j][col] != 0:
                    ma.MatRowSwap(M, j, row)
                    ind = True
                    break
            if ind == False :
                col += 1
                row += 1
    t = perf_counter() - t
    print("time execution: ",t)
    return M

#returns the inverse of a reversible matrix with Gauss-Jordan
def InverseMat(A):
    B = deepcopy(A)
    t = perf_counter()
    n = len(A)
    B = GaussPivot(B)
    B = ma.TransMat(B)
    B = GaussPivot(B)
    for c in range(n):
        B[c][c] = 1/B[c][c]
    t = perf_counter() - t
    print("time execution: ",t)
    return B
