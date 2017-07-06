from single_particle_state import SingleParticleState
class SingleParticleBasis(object):
    def __init__(this, states):
        this.states = states

    def __str__(this):
        s = ""
        for state in this.states:
            s += str(state) + "\n"
        return s[:-1]
