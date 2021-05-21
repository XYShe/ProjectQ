from ._basics import BasicGate,  SelfInverseGate
import numpy as np
from projectq.ops import get_inverse

class Reflection(BasicGate):

    def __init__(self, n):
        BasicGate.__init__(self)
        self.n = n

    def __str__(self):
        return 'Reflection Gate with n = {}'.format(str(self.n))

class SparseReflection(BasicGate):

    def __init__(self, vec,n):
        BasicGate.__init__(self)
        self.vec = vec
        self.n = n


    def __str__(self):
        return 'Sparse Householder Reflection({})'.format(str(self.vec))

class DenseHouseholderDec(BasicGate):

    def __init__(self, isometry):
        BasicGate.__init__(self)
        self.isometry = np.array(isometry,dtype='complex_')
    def __str__(self):
        return 'Dense Householder Decomposition({})'.format(str(self.isometry))



class SparseHouseholderDec(BasicGate):

    def __init__(self, isometry):
        BasicGate.__init__(self)
        self.isometry =  np.array(isometry,dtype='complex_')



    def __str__(self):
        return 'Sparse Householder Decomposition({})'.format(str(self.isometry))

