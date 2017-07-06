import numpy as np
from twobody_basis import TwobodyBasis

class TwobodyBlock(object):
    def __init__(this, basis, matrix=None):
        if type(basis) is TwobodyBasis:
            this.basis = basis
        else: this.basis = None

        this.shape = (this.basis.get_basis_size(), this.basis.get_basis_size())

        if type(matrix) is np.ndarray and len(matrix.shape) == 2:
            this.set_matrix(matrix)
        elif type(matrix) is list:
            this.set_matrix(np.ndarray(matrix))
        else:
            this.set_matrix(np.zeros(this.shape))

    def set_matrix(this, matrix):
        if this.shape == matrix.shape: this.matrix = matrix
        else: raise Exception("Matrix is not the correct size!")

    def get_matrix(this): return this.matrix

    def __str__(this):
        s = str(this.basis) + "\n"
        s += "Matrix shape: " + str(this.matrix.shape)
        return s

if __name__ == "__main__":
    from twobody_basis_creator import TwobodyBasisCreator
    from twobody_channel_creator import TwobodyChannelCreator
    from single_particle_basis_creator import SingleParticleBasisCreator as creator
    from itertools import compress

    spbasis = creator.create_sd_shell()
    twobody_basis = TwobodyBasisCreator(spbasis)
    channels = TwobodyChannelCreator.create_channels(spbasis)

    print "Proton proton channels"
    for ch in compress(channels, [ch.get_isospin() == -1 for ch in channels]):
        print TwobodyBlock(twobody_basis.get_basis_for_channel(ch))
        print

    print "Proton neutron channels"
    for ch in compress(channels, [ch.get_isospin() == 0 for ch in channels]):
        print TwobodyBlock(twobody_basis.get_basis_for_channel(ch))
        print

    print "Neutron neutron channels"
    for ch in compress(channels, [ch.get_isospin() == 1 for ch in channels]):
        print TwobodyBlock(twobody_basis.get_basis_for_channel(ch))
        print
