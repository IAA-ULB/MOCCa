"""
This file tests providing linearized and grid-based access to data.

see issues/51
"""
import numpy as np

class TWDA:
    def __init__(self, grid_shape, n_components):
        self.grid_shape = grid_shape
        self.linear_size = np.prod(grid_shape)
        self.n_components = n_components
        self.data = np.zeros((self.linear_size, self.n_components), dtype=int, order='F')
        self.grid_data = self.data.reshape((*self.grid_shape, n_components), order='F')

    def __getitem__(self, index):
        return self.data[index[1:]] if index[0] == 'L' else \
               self.grid_data[index[1:]]

def test_data_shared():
    """This test verifies that TWDA.data and TWDA.grid_data share the same data.
    Modifying one of them also modifies the other.
    """
    M = 2
    nq = 2
    twda = TWDA(grid_shape=(M,M), n_components=nq)

    # assign to data
    for iq in range(nq):
        for l in range(M*M):
            twda.data[l,iq] = l + 100*iq
    l = 0
    for iq in range(nq):
        for j in range(M):
            for i in range(M):
                print(f"({i},{j},{iq=}): G[{i},{j}] = {twda.grid_data[i,j,iq]} L({l=}) = {twda.data[l,iq]})")
                # grid_data and data should be equal
                assert twda.grid_data[i,j,iq] == twda.data[l,iq]
                l += 1
                l %= M*M

    # modify twda.grid_data and verify that twda.data changes accordingly
    # now the other way around:
    for iq in range(nq):
        for j in range(M):
            for i in range(M):
                twda.grid_data[i,j,iq] *= 2
    i,j = 0,0
    for iq in range(nq):
        for l in range(M*M):
            print(f"({i},{j},{iq=}): G[{i},{j}] = {twda.grid_data[i, j, iq]} L({l=}) = {twda.data[l, iq]})")
            assert twda.data[l,iq] % 2 == 0
            i += 1
            if i == M:
                i = 0
                j += 1
                if j == M:
                    j = 0

def test_getitem():
    M = 2
    nq = 2
    twda = TWDA(grid_shape=(M,M), n_components=nq)

    for iq in range(nq):
        for l in range(M*M):
            twda.data[l,iq] = l + 100*iq

    for iq in range(nq):
        for l in range(M*M):
            d = twda['L',l,iq]
            print(f"{iq=} L({l}) = {d}")
            assert d == l + 100*iq

    # slices work top
    print(twda['L', :, 0]) # print linear data for component 0
    print(twda['G', :, :, 0])  # print grid data for component 0