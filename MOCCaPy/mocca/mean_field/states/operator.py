import numpy as np

class Operator:
    def __init__(self, hfpsi):
        self.ket = hfpsi
        self.matrix = hfpsi.matrix()
        self.mesh = hfpsi.mesh

    def compute_derivatives(self):
        pass
    def local_factor(self):
        pass
    def scale(self):
        pass
    def matrix_representation(self):
        self.