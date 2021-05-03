from ._basics import BasicGate
import numpy as np

class Reflection(BasicGate):

    def __init__(self, n):
        BasicGate.__init__(self)
        self.n = n

    def __str__(self):
        return 'Reflection Gate with n = {}'.format(str(self.n))

class DenseHouseholderDec(BasicGate):

    def __init__(self, isometry):
        BasicGate.__init__(self)
        self.isometry = np.array(isometry,dtype='complex_')

    def __str__(self):
        return 'Dense Householder Decomposition({})'.format(str(self.isometry))


class SparseHouseholderDec(BasicGate):

    def __init__(self, isometry):
        BasicGate.__init__(self)
        self.isometry = isometry

    def __str__(self):
        return 'Sparse Householder Decomposition({})'.format(str(self.isometry))

