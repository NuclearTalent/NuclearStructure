from partial_wave import PartialWave
class SingleParticleState(object):
    isospin_projection = {
        -1 : "proton",
        1 : "neutron"
    }

    isospin_projection_reverse = {
        "proton" : -1,
        "neutron": 1
    }

    def __init__(this, n, l, j, tz):
        this.pw = PartialWave(n, l, j)
        this.isospin = SingleParticleState.isospin_projection[int(tz)]

    def __str__(this):
        return "%s: %s" % (this.get_isospin_symbol(), str(this.pw))

    def __eq__(this, other):
        return this.pw == other.pw and this.isospin == other.isospin

    def __neq__(this, other): return not this == other

    def get_isospin_symbol(this): return this.isospin
    def get_parity(this): return this.pw.get_parity()
    def get_isospin(this): return this.isospin_projection_reverse[this.isospin]
    def get_spin(this): return this.pw.get_spin()

if __name__ == "__main__":
    s = SingleParticleState(0, 0, 1, -1)
    print s
    s = SingleParticleState(3, 2, 3, 1)
    print s
