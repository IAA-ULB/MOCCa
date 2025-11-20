class GrammSchmidt:
    """Class for Gramm-Schmidt orthogonalisaton of HFPsi objects"""
    def __init__(self, A):
        """
        Args:
            hfpsi: wave function
            epsilon: di
        """
        self.hfpsi = hfpsi


    def orthogonalize(self, normalize=True):
        """Orthogonlize hfpsi.
        """
        overlap = Overlap(self.hfpsi)
        self.overlap_matrix = overlap.compute_matrix_representation()
        hfpsi_data3 = self.hfpsi.data.reshape((self.hfpsi.data.shape[0], 4, self.hfpsi.data.shape[1]//4), order='F')
        for ib in range(8):
            range_ib = self.hfpsi.hfblockrange[ib]
            if range_ib[1] > range_ib[0]:
                # block is not empty
                # Restrict all data structures to block ib
                epsilon_ib = self.epsilon[range_ib[0]:range_ib[1]]
                hfpsi_data3_ib = hfpsi_data3[:,:,range_ib[0]:range_ib[1]]
                overlap_ib = self.overlap_matrix[ib]
                order_ib = np.argsort(epsilon_ib)
                # for i in order_ib[1:]:
                #     for j in order_ib[0:i]:
                n = epsilon_ib.size
                for i in range(n):
                    for j in range(i):
                        hfpsi_data3[:,:,i] -= (overlap_ib[i,j] / overlap_ib[j,j]) * hfpsi_data3_ib[:,:,j]
        if normalize:
            self.normalize()

    def normalize(self):
        """Normalize hfpsi
        """
        overlap = Overlap(self.hfpsi)
        diagonal = overlap.compute_diagonal_elements()
        hfpsi_data3 = self.hfpsi.data.reshape((self.hfpsi.data.shape[0], 4, self.hfpsi.data.shape[1]//4), order='F')
        for i_spwf in range(self.hfpsi.n_total_wf):
            hfpsi_data3[:,:,i_spwf] /= diagonal[i_spwf]
