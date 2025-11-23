from pathlib import Path
import sys

def insert_path(rel_path='', verbose=False):
    """ Insert path MOCCaPy/rel_path in sys.path.
    E.g. to be able to import
      - tests (MOCCaPy/tests): call `insert_path()`
      - sandbox/*.py (MOCCaPy/sandbox/*.py): call `insert_path('sandbox')`
    """
    path = Path(__file__).parent
    while path.name != "MOCCaPy":
        path = path.parent
    path = path / rel_path
    if not path.exists():
        raise FileNotFoundError(f"Folder '{path}' does not exist.")
    sys.path.insert(0, str(path))
    if verbose:
        print(f"inserted path: {path}")