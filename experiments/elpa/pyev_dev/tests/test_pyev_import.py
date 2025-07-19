from pyev.helpers import *

def test_pyelpa_import():
    if use_pyelpa:
        from pyelpa import ProcessorLayout, DistributedMatrix, Elpa
    else:
        from pyev   import ProcessorLayout, DistributedMatrix, Elpa
    

def test_elpa():
    if use_pyelpa:
        from pyelpa import ProcessorLayout, DistributedMatrix, Elpa
    else:
        from pyev   import ProcessorLayout, DistributedMatrix, Elpa
    e = Elpa()

# ==============================================================================
# The code below is for debugging a particular test
# (normally all tests are run with pytest)
# ==============================================================================
if __name__ == "__main__":
    the_test_you_want_to_debug = test_pyelpa_import

    print("__main__ running", the_test_you_want_to_debug)
    the_test_you_want_to_debug()
    print("-*# finished #*-")
# ==============================================================================
