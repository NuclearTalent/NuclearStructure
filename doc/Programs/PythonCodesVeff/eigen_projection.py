import numpy as np

class EigenProjection(object):
    def project(this, pairs, mapping):
        return [p.project(mapping) for p in pairs]

if __name__ == "__main__":
    from eigensolver import EigenSolver, EigenSolverSymmetric
    from eigen_decomposition import EigenDecomposition
    from matrix_from_eigen_decomposition import MatrixFromEigenDecomposition

    a = np.identity(3)
    a[0,1] = 2; a[0,2] = 4; a[1,2] = 5
    a += a.transpose()

    solver = EigenSolverSymmetric()
    projection = EigenProjection()
    constructor = MatrixFromEigenDecomposition()
    decomposer = EigenDecomposition()

    pairs = solver.solve(a)

    mapping = np.array([0, 1])
    projected_pairs = projection.project(pairs[0:2], mapping)
    b = constructor.construct(*decomposer.decompose_from_pairs(projected_pairs))

    solver = EigenSolver()
    final_pairs = solver.solve(b)

    print "New eigenvalues"
    for p in final_pairs: print p.eigenvalue
    print
    print "Old eigenvalues"
    for p in pairs[0:2]: print p.eigenvalue


