import numpy as np

def matvec(A, x):
    m, n = A.shape
    y = np.zeros(x.shape)
    for i in range(m):
        for j in range(n):
            y[i] = y[i] + A[i,j]*x[j]
    return y
