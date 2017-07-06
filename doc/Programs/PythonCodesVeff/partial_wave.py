class PartialWave(object):
    spectral_notation = {
        "s" : 0,
        "p" : 1,
        "d" : 2,
        "f" : 3,
        "g" : 4,
        "h" : 5
    }

    spectral_notation_reverse = {
        0 : "s",
        1 : "p",
        2 : "d",
        3 : "f",
        4 : "g",
        5 : "h"
    }


    def __init__(this, n, l, j):
        this.n = int(n)
        this.l = int(l)
        this.j = int(j)

    def get_parity(this): return (-1)**this.l
    def get_spin(this): return this.j

    def __str__(this):
        return "%d%s_{%d/2}" % \
            ( this.n, this.spectral_notation_reverse[this.l], this.j)

    def __eq__(this, other):
        return this.n == other.n and this.l == other.l and this.j == other.j

    def __neq__(this, other): return not this == other

