import pycircle
import ctypes

# from pdb import set_trace as BREAKPOINT

print(f"\nstarting t.py")

circleso = pycircle.Wrapper('/data/antwerpen/201/vsc20170/tantalus_full/experiments/elpa/sandbox/libclass_circle.so')
circleso.class_circle_mp_test()
# circleso.class_circle_mp_test2()
print("\nfine?")

# BREAKPOINT()
p_circle = ctypes.c_void_p(circleso.alloc())
# print(f"{p_circle}")
# print(f"do we get here?")

circleso.set_radius(p_circle, ctypes.c_float(5.0))
circleso.print(p_circle)
# print(f"{ctypes.c_void_p(p_circle)=}")

# can we do sth with p_circle?

# circleso.class_circle_mp_set_radius(p_circle, ctypes.c_float(5.))
# circleso.class_circle_mp_circle_print(p_circle)

circleso.free(p_circle)