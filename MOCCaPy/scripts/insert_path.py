from pathlib import Path
import sys

def insert_path(rel_path=''):
    """ Insert path MOCCaPy/rel_path in sys.path."""
    path = Path(__file__).parent
    while path.name != "MOCCaPy":
        path = path.parent
    path = path / rel_path
    if not path.exists():
        raise FileNotFoundError(f"Folder '{path}' does not exist.")
    sys.path.insert(0, str(path))