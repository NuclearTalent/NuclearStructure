import numpy as np

class EigenPair(object):
    def __init__(this, eigenvalue, eigenvector):
        this.eigenvalue = eigenvalue
        this.eigenvector = np.array(eigenvector)

    def __str__(this):
        s = "Eigenvalue: %s\n" % this.eigenvalue
        s += "Eigenvector: %s" % str(this.eigenvector)
        return s

    def is_eigenpair(this, matrix):
        return np.allclose(np.matmul(matrix,this.eigenvector), this.eigenvalue*this.eigenvector)

    def project(this, mapping):
        eigenvalue = this.eigenvalue
        eigenvector = this.eigenvector[mapping]
        return EigenPair(eigenvalue, eigenvector)

