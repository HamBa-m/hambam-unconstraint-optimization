#returns sum of 2 matrices A and B
def MatSum(A , B):
    n = len(A)
    m = len(A[0])
    C = []
    for row in range(n):
        L=[]
        for col in range(m):
            L.append(A[row][col] + B[row][col])
        C.append(L)
    return C

#returns matrices product of A and B
def MatProd(A , B):
    n = len(A)
    m = len(B[0])
    p = len(A[0])
    C = []
    for row in range(n):
        L = []
        for col in range(m):
            sum = 0
            for k in range(p):
                sum += A[row][k] * B[k][col]
            L.append(sum)
        C.append(L)
    return C

#returns transposed matrice of A
def TransMat(A):
    n = len(A)
    m = len(A[0])
    B = []
    for i in range(m):
        l = []
        for j in range(n):
            l.append(A[j][i])
        B.append(l)
    return B

#swaps the rows i and j in M
def MatRowSwap(M ,i ,j):
    for k in range(len(M[0])):
        M[i][k], M[j][k] = M[j][k], M[i][k]
    return

#multiplies the row i of M by the coefficient coef
def MatRowDilatation(M, i, coef):
    for j in range(len(M[0])):
        M[i][j] *= coef
    return

#adds the coef*M[j] to the row M[i]
def MatRowTransvection(M, i, j, coef):
    for k in range(len(M[0])):
        M[i][k] += M[j][k] * coef
    return

#returns an (n*n) matrix of zeros
def ZerosMat(n):
    O = []
    for i in range(n):
        l = []
        for j in range(n):
            l.append(0)
        O.append(l)
    return O

#returns an (n*n) identity matrix
def IdentityMat(n):
    I = []
    for i in range(n):
        l = []
        for j in range(n):
            if i == j : l.append(1)
            else : l.append(0)
        I.append(l)
    return I

#prints a matrix
def MatPrint(M):
    for i in range(len(M)):
        print(M[i])
    print()
    return