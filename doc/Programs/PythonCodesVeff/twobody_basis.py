from twobody_channel import TwobodyChannel

class TwobodyBasis(object):
    def __init__(this, channel, states):
        if type(channel) is TwobodyChannel: this.channel = channel
        this.states = states

    def __str__(this):
        s = str(this.channel) +"\n"
        s += "Number of states: %d\n" % this.get_basis_size()
        for state in this.states:
            s += str(state) + "\n"
        return s[:-1]

    def get_basis_size(this): return len(this.states)
