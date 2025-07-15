import pycircle
import ctypes
circleso = pycircle.Wrapper('./libclass_circle.so')
circleso.class_circle_mp_test()
# circleso.class_circle_mp_test2()
print("fine?")

p = circleso.alloc()
print(f"do we get here? {p}")

circleso.free(p)