import numpy as np
import pytest

from mocca.mean_field.operators.overlap import Overlap


# @Timer("_projector()")
def _projector(hfpsi_d3_ib, i, j):
    """Used by gramm_schmidt."""
    Overlap_ij = np.einsum("hk,hk", hfpsi_d3_ib[:, :, i], hfpsi_d3_ib[:, :, j], order='F', optimize=True)
    Overlap_jj = np.einsum("hk,hk", hfpsi_d3_ib[:, :, j], hfpsi_d3_ib[:, :, j], order='F', optimize=True)
    # both are missing a factor mesh.dv but that doesn's matter because of the quotient
    return Overlap_ij / Overlap_jj

# @Timer('gramm_schmidt()')
def gramm_schmidt(hfpsi, order=None, normalize=True, check=False):
    """Gramm-Schmidt orthogonalization of a wave function.

    Args:
        hfpsi: wave function to orthogonalize
        order: ndarray of length hfpsi.n_total_wf. This array determines the selection order of the
            single particle wave functions uin the orthogonalization process. The lowest spwf is
            selected first. If None, the order is the order of the single particle wave functions
            in hfpsi.data.
        normalize: whether to normalize the single particle wave functions after orthogonalization.
        check: if True, verifies that the overlap matrix is unity. (For debugging purposes only).
    """
    hfpsi_d3 = hfpsi.d3
    for ib in range(8):
        range_ib = hfpsi.hfblockrange[ib]
        if range_ib[1] > range_ib[0]:
            # block is not empty
            # Restrict all data structures to symmetry block ib
            hfpsi_d3_ib = hfpsi_d3[:,:,range_ib[0]:range_ib[1]]
            if not order is None:
                order_ib = np.argsort(order[range_ib[0]:range_ib[1]])
            n = range_ib[1]-range_ib[0]
            for I in range(1,n):
                i = order_ib[I] if (order is not None) else I
                for J in range(I):
                    j = order_ib[J] if (order is not None) else J
                    hfpsi_d3_ib[:,:,i] -= _projector(hfpsi_d3_ib,i,j) * hfpsi_d3_ib[:,:,j]
                # for J in range(0,I):
                #     j = order_ib[J]
                #     Oij = np.einsum("hk,hk", hfpsi_d33_ib[:,:,i], hfpsi_d33_ib[:,:,j])
                #     print(f"{ib=} ({i},{j}) {Oij=}")

    if normalize:
        hfpsi.normalize()

    if check:
        overlap = Overlap(hfpsi)
        overlap.compute_matrix_representation()
        for ib in range(8):
            n_ib = hfpsi.hfblocks[ib]
            block_ib = overlap.operand_data.matrix[ib]
            for i in range(0, n_ib):
                for j in range(i, n_ib):
                    if i == j:
                        if normalize:
                            assert block_ib[i,i] == pytest.approx(1.0), f"block[{ib}][{i},{j}] {block_ib[i, j]}"
                    else:
                        assert block_ib[i,j] == pytest.approx(0.0), f"block[{ib}][{i},{j}] {block_ib[i, j]}"