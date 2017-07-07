import numpy as np

class MatrixFromEigenDecomposition(object):
    def construct(this, D, V):
        return np.matmul(V, np.matmul(D, np.linalg.inv(V)))

class MatrixFromEigenDecompositionSymmetric(object):
    def construct(this, D, V):
        return np.matmul(V, np.matmul(D, V.transpose()))
