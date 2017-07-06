from single_particle_state import SingleParticleState
from single_particle_basis import SingleParticleBasis
class SingleParticleBasisCreator(object):
    pshell_partial_waves = [ 
        [ 0, 0, 1 ],
        [ 0, 1, 3 ],
        [ 0, 1, 1 ] ]

    sdshell_partial_waves = [
        [ 0, 2, 5 ],
        [ 1, 0, 1 ],
        [ 0, 2, 3 ] ]

    @staticmethod
    def create_p_shell(isospin=None):
        return SingleParticleBasisCreator.create_shell(
            SingleParticleBasisCreator.pshell_partial_waves, isospin)

    @staticmethod
    def create_sd_shell(isospin=None):
        return SingleParticleBasisCreator.create_shell(
            SingleParticleBasisCreator.sdshell_partial_waves, isospin)

    @staticmethod
    def create_shell(partial_waves, isospin):
        if isospin is None: isospin = [-1, 1 ]
        else: isospin = [isospin]

        states = []
        for tz in isospin:
            for pw in partial_waves:
                n, l, j = pw
                states.append(SingleParticleState(n, l, j, tz))
        return SingleParticleBasis(states)


if __name__ == "__main__":
    print "sd-shell states"
    sd = SingleParticleBasisCreator.create_sd_shell()
    print sd

    print
    print "p-shell states"
    p = SingleParticleBasisCreator.create_p_shell()
    print p
