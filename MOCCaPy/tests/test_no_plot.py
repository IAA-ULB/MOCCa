
import inspect

import pytest
from mocca.util import started_finished, started, finished


def test_answer(no_plot, debug=False):
    if no_plot:
        print("\nnot plotting\n")
    else:
        print("\nplotting\n")

# This one FAILS. The `started_finished` decorator causes the fixture noplot not to resolve.
@started_finished
def answer_decorated(no_plot):
    if no_plot:
        print("\nnot plotting\n")
    else:
        print("\nplotting\n")

def test_answer_decorated():
    with pytest.raises(TypeError):
        answer_decorated()

def test_answer_revisited(no_plot, debug=False):
    started(inspect.stack()[0][3])

    if no_plot:
        print("\nnot plotting\n")
    else:
        print("\nplotting\n")

    finished(inspect.stack()[0][3])