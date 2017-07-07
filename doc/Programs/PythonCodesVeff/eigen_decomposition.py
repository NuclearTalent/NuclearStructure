from eigensolver import EigenSolver, EigenSolverSymmetric
import numpy as np
class EigenDecomposition(object):
    def __init__(this):
        this.solver = EigenSolver()

    def decompose(this, matrix):
        return this.decompose_from_pairs(this.solver.solve(matrix))

    def decompose_from_pairs(this, pairs):
        shape = (len(pairs), len(pairs))
        D = np.zeros(shape)
        V = np.zeros(shape)
        D[np.diag_indices_from(D)] = np.array([p.eigenvalue for p in pairs])
        for i in range(len(pairs)): V[:,i] = pairs[i].eigenvector
        
        return D, V

class EigenDecompositionSymmetric(EigenDecomposition):
    def __init__(this):
        this.solver = EigenSolverSymmetric()

if __name__ == "__main__":
    from matrix_from_eigen_decomposition import MatrixFromEigenDecomposition, \
        MatrixFromEigenDecompositionSymmetric

    constructor = MatrixFromEigenDecomposition()
    decomposer = EigenDecomposition()

    a = np.identity(3)
    a[0,1] = 2; a[0,2] = 4; a[1,2] = 5
    a += a.transpose()
    
    D, V = decomposer.decompose(a)
    ap = constructor.construct(D, V)
    if not np.allclose(a, ap):
        res = "FAILED"
    else: res = "PASSED"
    print "Test EigenDecomposition; %s" % res

    decomposer = EigenDecompositionSymmetric()
    constructor = MatrixFromEigenDecompositionSymmetric()
    D, V = decomposer.decompose(a)
    ap = constructor.construct(D, V)
    if not np.allclose(a, ap):
        res = "FAILED"
    else: res = "PASSED"
    print "Test EigenDecompositionSymmetric; %s" % res
