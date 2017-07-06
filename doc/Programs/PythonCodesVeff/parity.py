class Parity(object):
    def __init__(this, parity):
        if type(parity) is int: this.parity = parity
        elif type(parity) is str: 
            if parity == "+": this.parity = 1
            elif parity == "-": this.parity = -1
        elif type(parity) is Parity: this.parity = parity.get_parity()

    def get_parity_symbol(this):
        if this.parity == 1: return "+"
        elif this.parity == -1: return "-"

    def get_parity(this): return this.parity
