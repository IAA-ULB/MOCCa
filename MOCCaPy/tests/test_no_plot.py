
from mocca.util import started_finished

def test_answer(no_plot, debug=False):
    if no_plot:
        print("\nnot plotting\n")
    else:
        print("\nplotting\n")

# This one FAILS. The decorator causes the fixture noplot not to resolve.
# @started_finished
# def test_answer_decorated(no_plot, debug=False):
#     if no_plot:
#         print("\nnot plotting\n")
#     else:
#         print("\nplotting\n")
