import math
import numpy as np

from projectq.cengines import DecompositionRule
from projectq.meta import Control, Compute, Uncompute, CustomUncompute, Dagger
from projectq.ops import X, Z, Ph, All, StatePreparation, StatePrepRecur, \
    DaggeredGate, Rx,Ry,Rz, CNOT, C, MatrixGate, SparseReflection
from projectq.libs._utils import schmidtDec, ComputeHouseholderVec, chop, stdHHref, \
    isPermutedDiag, decPermDiag, applyGatetoISO, nextGen
from projectq.ops import DenseHouseholderDec, SparseHouseholderDec, Reflection
from scipy import sparse





def _decompose_reflect(cmd):
    eng = cmd.engine
    n = cmd.gate.n
    qureg = cmd.qubits[0]
    with Control(eng, cmd.control_qubits):
        Rz(np.pi) | qureg[n-1]
        Ry(np.pi/2) | qureg[n-1]
        Rz(np.pi) | qureg[n-1]

        for i in range(n-1):
            Rz(np.pi) | qureg[i]
            Ry(np.pi) | qureg[i]

        if n == 1:
            Rz(np.pi) | qureg[0]
            Ry(np.pi) | qureg[0]


        elif n == 2:
            CNOT | (qureg[0], qureg[1])

        elif n >= 3:
            print('Invoked')
            ctrl_list = []

            for i in range(n-1):
                ctrl_list.append(qureg[i])
            C(X,n-1) | (ctrl_list, qureg[n])

        for i in range(n-1):
            Rz(np.pi) | qureg[i]
            Ry(np.pi) | qureg[i]

        #########################3
        Ry(np.pi/2) | qureg[n-1]




def _decompose_sparse_reflect(cmd):
    eng = cmd.engine
    n = cmd.gate.n
    qureg = cmd.qubits[0]
    vec = cmd.gate.vec
    with Control(eng, cmd.control_qubits):
        DaggeredGate(StatePreparation(vec)) | qureg
        Reflection(n) | qureg
        StatePreparation(vec) | qureg



def _decompose_dense(cmd):
    eng = cmd.engine
    assert len(cmd.qubits) == 1
    num_qubits = len(cmd.qubits[0])
    qureg = cmd.qubits[0]

    V = cmd.gate.isometry
    assert np.log2(V.shape[1]).is_integer() == True, 'Dimension should be power of 2'
    assert np.log2(V.shape[0]).is_integer() == True, 'Dimension should be power of 2'

    m = int(np.log2(V.shape[1]))
    n = int(np.log2(V.shape[0]))


    assert m == num_qubits, 'The number of columns should equal to number of qubits'
    assert m <= n, 'Isometry requires dim1 <= dim2'

    with Control(eng, cmd.control_qubits):

        if m < n:
            aux = eng.allocate_qureg(n-m)
            qureg.extend(aux)

        vecs = []

        for i in range(2**m):
            v = np.array(V[:,[i]])

            vec = ComputeHouseholderVec(v, i)

            if np.abs(vec[i]) != 1:

                V = stdHHref(vec,V)
                vecs.append(vec)

        #diagonal =np.diag(np.diag(V))
        #MatrixGate(diagonal) | qureg
        print(vecs[-1])
        DaggeredGate(StatePreparation(vecs[-1]))| qureg



        for i in range(1,len(vecs)):
            print('invoked')
            Reflection(n) | qureg
            DaggeredGate(StatePreparation(vecs[i]))|qureg
            StatePreparation(vecs[i-1]) | qureg


        Reflection(n) | qureg
        StatePreparation(vecs[0])| qureg

        #diagonal =np.diag(np.diag(V))
        #MatrixGate(diagonal) | qureg
        #print(diagonal)

def _decompose_sparse(cmd):
    eng = cmd.engine
    assert len(cmd.qubits) == 1
    num_qubits = len(cmd.qubits[0])
    qureg = cmd.qubits[0]

    iso = cmd.gate.isometry
    print("Sparsed")

    assert np.log2(iso.shape[1]).is_integer() == True, 'Dimension should be power of 2'
    assert np.log2(iso.shape[0]).is_integer() == True, 'Dimension should be power of 2'

    m = int(np.log2(iso.shape[1]))
    n = int(np.log2(iso.shape[0]))


    assert m == num_qubits, 'The number of columns should equal to number of qubits'
    assert m <= n, 'Isometry requires dim1 <= dim2'


    with Control(eng, cmd.control_qubits):
        if m < n:
            aux = eng.allocate_qureg(n-m)
            qureg.extend(aux)

        V = sparse.csr_matrix(iso)
        gatelist = []
        counter = 0
        while not isPermutedDiag(V):
            (i,j) = nextGen(V)
            vec = ComputeHouseholderVec(np.copy(iso[:,j]),i)
            SparseReflection(vec,num_qubits) | qureg
            vec = np.reshape(vec,(1,vec.shape[0]))
            iso = stdHHref(vec.transpose(),iso)
            V = sparse.csr_matrix(iso)
            counter += 1



        if isPermutedDiag(V):
            gates = decPermDiag(V)
            iso = applyGatetoISO(gates, iso)
            gatelist.extend(gates)
            diagonal =np.diag(np.diag(iso))
            gatelist.append(diagonal)

        for gate in gatelist[::-1]:
            MatrixGate(gate) | qureg

#: Decomposition rules
all_defined_decomposition_rules = [DecompositionRule(Reflection, _decompose_reflect),
DecompositionRule(SparseReflection, _decompose_sparse_reflect),
DecompositionRule(DenseHouseholderDec, _decompose_dense),
DecompositionRule(SparseHouseholderDec, _decompose_sparse)]


