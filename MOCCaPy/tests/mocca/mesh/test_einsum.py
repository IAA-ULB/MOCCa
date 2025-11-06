import numpy as np

def test_einsum_1D():
    D = np.array([[1, 2],
                  [3, 4]], dtype=float, order='F')
    # 1 component
    q = np.array([5, 6], dtype=float, order='F')
    Dq = np.einsum('il,l', D, q)
    assert Dq.shape == (2, )
    assert Dq[0] == 1*5 + 2*6
    assert Dq[1] == 3*5 + 4*6
    print(f"{Dq=}")

    # 2 components
    q = np.array([[5,5],
                  [6,6]], dtype=float, order='F')
    Dq = np.einsum('il,lq', D, q)
    assert Dq.shape == (2, 2)
    for i in range(2):
        assert Dq[0,i] == 1*5 + 2*6
        assert Dq[1,i] == 3*5 + 4*6
    print(f"{Dq=}")

    # 3 components
    q = np.array([[5,5,5],
                  [6,6,6]], dtype=float, order='F')
    Dq = np.einsum('il,lq', D, q)
    assert Dq.shape == (2, 3)
    for i in range(3):
        assert Dq[0,i] == 1*5 + 2*6
        assert Dq[1,i] == 3*5 + 4*6
    print(f"{Dq=}")
    
def test_einsum_2D():
    D = np.array([[1, 2],
                  [3, 4]], dtype=float, order='F')
    print(f"{D=}")
    # 1 component
    q = np.array([[1, ], # (0,0)
                  [2, ], # (1,0)
                  [3, ], # (0,1)
                  [4, ], # (1,1)
                  ], dtype=float, order='F').reshape((2,2), order='F') # reshaping to grid dimensions,
    print(f"{q=}")
    Dqx = np.einsum('il,lj', D, q)
    Dqy = np.einsum('il,jl', D, q)
    print(f"{Dqx=}")
    print(f"{Dqy=}")
    assert Dqx.shape == (2, 2)
    assert Dqy.shape == (2, 2)

    assert Dqx[0,0] == 1*1 + 2*2
    assert Dqx[1,0] == 3*1 + 4*2
    assert Dqx[0,1] == 1*3 + 2*4
    assert Dqx[1,1] == 3*3 + 4*4

    assert Dqy[0,0] == 1*1 + 2*3
    assert Dqy[1,0] == 3*1 + 4*3
    assert Dqy[0,1] == 1*2 + 2*4
    assert Dqy[1,1] == 3*2 + 4*4

    # 2 componentS
    q = np.array([[1, 1],  # (0,0)
                  [2, 2],  # (1,0)
                  [3, 3],  # (0,1)
                  [4, 4],  # (1,1)
                  ], dtype=float, order='F').reshape((2, 2, 2), order='F') # reshaping to grid/q dimensions
    print(f"{q=}")
    Dqx = np.einsum('il,ljq', D, q)
    Dqy = np.einsum('il,jlq', D, q)
    print(f"{Dqx=}")
    print(f"{Dqy=}")
    assert Dqx.shape == (2, 2, 2)
    assert Dqy.shape == (2, 2, 2)
    for i in range(2):
        assert Dqx[0, 0, i] == 1 * 1 + 2 * 2
        assert Dqx[1, 0, i] == 3 * 1 + 4 * 2
        assert Dqx[0, 1, i] == 1 * 3 + 2 * 4
        assert Dqx[1, 1, i] == 3 * 3 + 4 * 4
        assert Dqy[0, 0, i] == 1 * 1 + 2 * 3
        assert Dqy[1, 0, i] == 3 * 1 + 4 * 3
        assert Dqy[0, 1, i] == 1 * 2 + 2 * 4
        assert Dqy[1, 1, i] == 3 * 2 + 4 * 4

    