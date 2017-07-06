from single_particle_state import SingleParticleState
from parity import Parity

class TwobodyState(object):
    def __init__(this, a, b):
        if type(a) is SingleParticleState and type(b) is SingleParticleState:
            this.a = a
            this.b = b
        else: this.a = None; this.b = None

    def is_identical(this): return this.a == this.b

    def __str__(this):
        return "a: %s\nb: %s\nIdentical: %s" %(str(this.a), str(this.b),this.is_identical())

    def get_parity_symbol(this): return Parity(this.a.get_parity()*this.b.get_parity()).get_parity_symbol()
    def get_parity(this): return Parity(this.a.get_parity()*this.b.get_parity()).get_parity()
    def get_isospin(this): return (this.a.get_isospin() + this.b.get_isospin())/2
    def get_minimum_spin(this): return abs(this.a.get_spin() - this.b.get_spin())/2
    def get_maximum_spin(this): return (this.a.get_spin() + this.b.get_spin())/2

if __name__ == "__main__":
    t = TwobodyState(SingleParticleState(0, 2, 5, 1), SingleParticleState(0, 2, 5, 1))
    print t
    t = TwobodyState(SingleParticleState(0, 2, 5, 1), SingleParticleState(0, 2, 5, -1))
    print t
