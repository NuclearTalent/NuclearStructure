from twobody_channel import TwobodyChannel
from twobody_state import TwobodyState
from twobody_basis import TwobodyBasis

class TwobodyBasisCreator(object):
    def __init__(this, spbasis):
        this.spbasis = spbasis

    def get_basis_for_channel(this, channel):
        if type(channel) is not TwobodyChannel: return []

        states = []
        left = this.spbasis.states
        right = this.spbasis.states
        for a, idxa in zip(left, range(len(left))):
            for b, idxb in zip(right, range(len(right))):
                if idxb < idxa: continue
                state = TwobodyState(a, b)
                if state.get_parity() != channel.get_parity(): continue
                if state.get_isospin() != channel.get_isospin(): continue
                if channel.get_spin() < state.get_minimum_spin(): continue
                if channel.get_spin() > state.get_maximum_spin(): continue
                if a==b and channel.get_spin()%2 != 0: continue
                states.append(state)
        return TwobodyBasis(channel, states)

if __name__ == "__main__":
    from twobody_channel_creator import TwobodyChannelCreator
    from single_particle_basis_creator import SingleParticleBasisCreator as creator
    from itertools import compress

    spbasis = creator.create_sd_shell()
    twobody_basis = TwobodyBasisCreator(spbasis)
    channels = TwobodyChannelCreator.create_channels(spbasis)

    print "Proton proton channels"
    for ch in compress(channels, [ch.get_isospin() == -1 for ch in channels]):
        print twobody_basis.get_basis_for_channel(ch)
        print
    
    print "Proton neutron channels"
    for ch in compress(channels, [ch.get_isospin() == 0 for ch in channels]):
        print twobody_basis.get_basis_for_channel(ch)
        print
    
    print "Neutron neutron channels"
    for ch in compress(channels, [ch.get_isospin() == 1 for ch in channels]):
        print twobody_basis.get_basis_for_channel(ch)
        print

