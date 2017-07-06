import numpy as np
from eigenpair import EigenPair
class EigenSolver(object):
    def solve(this, matrix):
        w, v = np.linalg.eig(matrix)
        return [EigenPair(w[i], v[:,i]) for i in range(len(w))]

class EigenSolverSymmetric(EigenSolver):
    def solve(this, matrix):
        w, v = np.linalg.eigh(matrix)
        return [EigenPair(w[i], v[:,i]) for i in range(len(w))]

if __name__ == "__main__":
    solver = EigenSolver()

    a = np.identity(3)
    a[0,1] = 2; a[0,2] = 4; a[1,2] = 5
    a += a.transpose()

    pairs = solver.solve(a)
    
    if not all([ p.is_eigenpair(a) for p in pairs]):
        res = "FAILED"
    else: res = "PASSED"
    print "Test EigenSolver; %s" % res

    solver = EigenSolver()

    a = np.identity(3)
    a[0,1] = 2; a[0,2] = 4; a[1,2] = 5
    a += a.transpose()

    pairs = solver.solve(a)
    if not all([ p.is_eigenpair(a) for p in pairs]):
        res = "FAILED"
    else: res = "PASSED"
    print "Test EigenSolverSymmetric; %s" % res
