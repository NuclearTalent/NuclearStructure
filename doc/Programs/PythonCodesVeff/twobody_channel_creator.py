from twobody_channel import TwobodyChannel
from parity import Parity

class TwobodyChannelCreator(object):
    @staticmethod
    def create_channels(basis):
        channels = []
        parity_channels = TwobodyChannelCreator.find_parity(basis)
        isospin_channels = TwobodyChannelCreator.find_isospin(basis)
        spin_channels = TwobodyChannelCreator.find_spin(basis)
        for parity in parity_channels:
            for isospin in isospin_channels:
                for spin in spin_channels:
                    channels.append(TwobodyChannel(spin, isospin, parity))

        return channels

    @staticmethod
    def find_parity(basis):
        parity = [s.get_parity() for s in basis.states]
        if abs(sum(parity)) == len(parity): return [ Parity('+') ]
        else: return [Parity('+'), Parity("-")]

    @staticmethod
    def find_isospin(basis):
        isospin = [s.get_isospin() for s in basis.states]
        if sum(isospin) == len(isospin): return [1]
        elif abs(sum(isospin)) == len(isospin): return [-1]
        else: return [-1, 0, 1]

    @staticmethod
    def find_spin(basis): return range(max([s.get_spin() for s in basis.states]) + 1)

if __name__ == "__main__":
    from single_particle_basis_creator import SingleParticleBasisCreator

    print "p-shell"
    p = SingleParticleBasisCreator.create_p_shell()
    channels = TwobodyChannelCreator.create_channels(p)
    for c in channels: print c

    print "sd-shell"
    sd = SingleParticleBasisCreator.create_sd_shell()
    channels = TwobodyChannelCreator.create_channels(sd)
    for c in channels: print c
