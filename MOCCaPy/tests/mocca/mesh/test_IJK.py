from mocca.mesh.ijk import IJK
import pytest


def test_IJK_1():
    print("test_IJK_1")
    shape = (3,)
    ijk = IJK(shape)
    i,_,_ = ijk.index
    for ig in range(shape[0]):
        print(f"{ig=} : {ijk.index}")
        assert ijk.index[0] == ig
        assert i == ig
        ijk_from_ig = ijk.ig2ijk(ig)
        assert ijk_from_ig == ijk.index
        ig_from_ijk = ijk.ijk2ig(ijk_from_ig)
        assert ig_from_ijk == ig

        i,_,_ = ijk.inc()


def test_IJK_2():
    print("test_IJK_2")
    shape = (3,2)
    ijk = IJK(shape)
    i,j,_ = ijk.index
    for ig in range(shape[0]*shape[1]):
        print(f"{ig=} : {ijk.index}")
        assert ijk.index[0] +ijk.index[1]*shape[0] == ig
        assert i + j*shape[0] == ig
        ijk_from_ig = ijk.ig2ijk(ig)
        assert ijk_from_ig == ijk.index
        ig_from_ijk = ijk.ijk2ig(ijk_from_ig)
        assert ig_from_ijk == ig

        i,j,_ = ijk.inc()

def test_IJK_2_done():
    print("test_IJK_2_done")
    shape = (3,2)
    ijk = IJK(shape)
    while not ijk.done():
        print(f"{ijk.ig} : {ijk.index}")
        i,j,_ = ijk.inc()

def test_IJK_3():
    print("test_IJK_3")
    shape = (4, 3, 2)
    ijk = IJK(shape)
    i,j,k = ijk.index
    for ig in range(shape[0] * shape[1] * shape[2]):
        print(f"{ig=} : {ijk.index}")
        assert ijk.index[0] + ijk.index[1]*shape[0] + ijk.index[2]*(shape[0] * shape[1]) == ig
        assert i + j*shape[0] + k*(shape[0] * shape[1]) == ig
        ijk_from_ig = ijk.ig2ijk(ig)
        assert ijk_from_ig == ijk.index
        ig_from_ijk = ijk.ijk2ig(ijk_from_ig)
        assert ig_from_ijk == ig

        i,j,k = ijk.inc()

def test_IJK_3_done():
    print("test_IJK_3_done")
    shape = (4,3,2)
    ijk = IJK(shape)
    while not ijk.done():
        print(f"{ijk.ig} : {ijk.index}")
        i,j,k = ijk.inc()
