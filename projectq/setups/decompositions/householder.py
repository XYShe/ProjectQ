import math
import numpy as np

from projectq.cengines import DecompositionRule
from projectq.meta import Control, Compute, Uncompute, CustomUncompute, Dagger
from projectq.ops import X, Z, Ph, All, StatePreparation, DaggeredGate, Rx,Ry,Rz, CNOT, C, MatrixGate
from projectq.libs._utils import SchmidtDec, ComputeHouseholderVec
from projectq.ops import DenseHouseholderDec, SparseHouseholderDec, Reflection





def _decompose_reflect(cmd):
    eng = cmd.engine
    n = cmd.gate.n
    qureg = cmd.qubits[0]
    with Control(eng, cmd.control_qubits):
        Ry(np.pi/2) | qureg[n-1]
        for i in range(n-1):
            Ry(np.pi) | qureg[i]
            Rz(np.pi) | qureg[i]

        if n == 1:
            Ry(np.pi) | qureg[0]
            Rz(np.pi) | qureg[0]

        elif n == 2:
            CNOT | (qureg[0], qureg[1])

        elif n >= 3:
            print('Invoked')
            ctrl_list = []

            for i in range(n-1):
                ctrl_list.append(qureg[i])
            C(X,n-1) | (ctrl_list, qureg[n])

        for i in range(n-1):
            Ry(np.pi) | qureg[i]
            Rz(np.pi) | qureg[i]

        Rz(np.pi) | qureg[n-1]
        Ry(np.pi/2) | qureg[n-1]
        Rz(np.pi) | qureg[n-1]

def _decompose_dense(cmd):
    eng = cmd.engine
    assert len(cmd.qubits) == 1
    num_qubits = len(cmd.qubits[0])
    qureg = cmd.qubits[0]

    V = cmd.gate.isometry
    print(V.shape)
    assert np.log2(V.shape[1]).is_integer() == True
    assert np.log2(V.shape[0]).is_integer() == True
    assert np.log2(V.shape[1]) == num_qubits

    with Control(eng, cmd.control_qubits):
        m = int(np.log2(V.shape[1]))
        n = int(np.log2(V.shape[0]))
        vecs = []

        for i in range(2**m):

            vec = ComputeHouseholderVec(V[:,i], i)
            vecs.append(vec)
        print(vecs)
        DaggeredGate(StatePreparation(vecs[0])) | qureg
        Reflection(n) | qureg

        for i in range(1,len(vecs)):
            StatePreparation(vecs[i-1]) | qureg
            DaggeredGate(StatePreparation(vecs[i]))|qureg
            Reflection(n) | qureg

        DaggeredGate(StatePreparation(vecs[-1])) | qureg

        diagonal =np.diag(np.diag(V))
        MatrixGate(diagonal) | qureg

def _decompose_sparse(cmd):
    """ Decompose the Quantum Amplitude Apmplification algorithm as a gate. """
    eng = cmd.engine

    # System-qubit is the first qubit/qureg. Ancilla qubit is the second qubit
    #system_qubits = cmd.qubits[0]
    pass

#: Decomposition rules
all_defined_decomposition_rules = [DecompositionRule(Reflection, _decompose_reflect),
DecompositionRule(DenseHouseholderDec, _decompose_dense),
DecompositionRule(SparseHouseholderDec, _decompose_sparse)]


