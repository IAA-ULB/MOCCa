from .util.timer import Timer
from pathlib import Path
import json

__version__ = "3.0.0"

package_folder = Path(__file__)
data_folder = Path(__file__).parent / "data"
data_folder.mkdir(parents=True, exist_ok=True)
if (data_folder/'chemical_symbols.json').exists():
    with open(data_folder/'chemical_symbols.json', 'r') as f:
        chemical_symbols = json.load(f)
else:
    # `mendeleev` package provides lots of data on elements/isotopes. I wanted a conversion from the atomic
    # number to the element's symbol. `mendeleev` does have this (as well as a lot of other stuff, it might
    # be perceived as overkill..., also because it has a lot of dependencies that we do not really need.)
    # see https://mendeleev.readthedocs.io/en/stable/index.html.
    try:
        from mendeleev.fetch import fetch_table
        elements = fetch_table('elements')
        chemical_symbols = [None, ]
        try:
            for Z in range(1,200):
                symbol = elements['symbol'][Z - 1]
                # print(Z, symbol)
                chemical_symbols.append(symbol)
        except KeyError:
            pass
        with open(data_folder /'chemical_elements.json', 'w') as f:
            json.dump(chemical_symbols, f)
    except ModuleNotFoundError:
        pass


def symbol(Z: int) -> str:
    """Return the chemical symbol corresponding to the given atomic number Z."""
    # The row index of the table starts at 0, hence the row index needed is Z-1.
    return chemical_symbols[Z]


