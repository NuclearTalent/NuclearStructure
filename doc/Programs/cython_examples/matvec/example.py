from matvec import matvec
import numpy as np
import cProfile as p

def run_cython(A, x):
    return matvec(A, x)

def run_python(A, x):
    m, n = A.shape
    y = np.zeros(x.shape)
    for i in range(m):
        for j in range(n):
            y[i] = y[i] + A[i,j]*x[j]

    return y

def run_blas(A, x):
    return np.matmul(A,x)

def run_vectorized(A, x):
    return np.array([(A[i, :]*x).sum() for i in range(A.shape[1])])

if __name__ == "__main__":
    import sys

    usage = "Usage: python %s dim" % sys.argv[0]
    try:
        m = int(sys.argv[1])
    except:
        print usage
        sys.exit(1)

    A = np.random.random_sample((m,m))
    x = np.ones((m))

    result = np.allclose(run_blas(A, x), run_cython(A, x)) and \
        np.allclose(run_blas(A, x), run_python(A, x)) and \
        np.allclose(run_blas(A, x), run_vectorized(A, x))
    if not result: print "matvec functions does not match!"

    p.run("run_python(A, x)")
    p.run("run_cython(A, x)")
    p.run("run_vectorized(A, x)")
    p.run("run_blas(A, x)")
