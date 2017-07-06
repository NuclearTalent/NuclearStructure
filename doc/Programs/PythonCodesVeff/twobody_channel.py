from parity import Parity
class TwobodyChannel(object):
    isospin_conversion = {
        -1 : "pp",
        0 : "pn",
        1 : "nn"
    }

    def __init__(this, spin, isospin, parity):
        this.spin = int(spin)
        this.isospin = int(isospin)
        this.parity = Parity(parity)

    def get_spin(this): return this.spin
    def get_isospin(this): return this.isospin
    def get_parity(this): return this.parity.get_parity()

    def __str__(this):
        return "Spin: %d, Isospin: %s, Parity: %s" % (this.spin, this.get_isospin_symbol(), this.get_parity_symbol())

    def get_parity_symbol(this): return this.parity.get_parity_symbol()

    def get_isospin_symbol(this):
        return this.isospin_conversion[this.isospin]

if __name__ == "__main__":
    c = TwobodyChannel(0, 0, '+')
    print c, c.get_spin(), c.get_isospin(), c.get_parity()
    c = TwobodyChannel(3, 1, '-')
    print c, c.get_spin(), c.get_isospin(), c.get_parity()
    c = TwobodyChannel(4, -1, '+')
    print c, c.get_spin(), c.get_isospin(), c.get_parity()
