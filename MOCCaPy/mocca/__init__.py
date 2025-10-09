__version__ = "3.0.0"

from mendeleev.fetch import fetch_table
elements = fetch_table('elements')

def symbol(Z: int) -> str:
    """Return the chemical symbol corresponding to the given atomic number Z."""
    # The row index of the table starts at 0, hence the row index needed is Z-1.
    return elements['symbol'][Z-1]

class Nucleus:
    def __init__(self, Z: int, N: int):
        self.Z = Z # number of protons (atomic number)
        self.N = N # number of neutrons
        self.symbol = symbol(Z)

    def __repr__(self):
        return f'Nucleus({self.symbol}, Z={self.Z}, N={self.N})'

    def __str__(self):
        return f'{self.symbol}[{self.Z},{self.N}]'