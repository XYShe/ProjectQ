import pytest
from cmath import exp
import math

import numpy as np
import projectq
from projectq.backends import Simulator
from projectq.cengines import (AutoReplacer, DecompositionRuleSet,
                               DummyEngine, InstructionFilter, MainEngine)
from projectq.ops import (BasicGate, ClassicalInstructionGate, MatrixGate,
                          Measure, Ph, R, Rx, Ry, Rz, X, DenseHouseholderDec, SparseHouseholderDec, All)
from projectq.meta import Control

from scipy.stats import unitary_group
from projectq.setups.decompositions import householder
import projectq.libs.math
from projectq.libs._utils import create_iso
def create_test_matrices():
    #matrix = unitary_group.rvs(4)
    return [np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])]




@pytest.mark.parametrize("gate_matrix", create_iso(2,3))
def test_dense_decomposition(gate_matrix):
    for basis_state in ([1, 0,0,0], [0,0,0, 1]):
        # Create single qubit gate with gate_matrix
        n=2
        test_gate = MatrixGate()
        test_gate.matrix = np.matrix(gate_matrix)

        correct_dummy_eng = DummyEngine(save_commands=True)
        correct_eng = MainEngine(backend=Simulator(),
                                 )
        rule_set = DecompositionRuleSet(modules=[projectq.libs.math,
                                                 projectq.setups.decompositions])
        test_dummy_eng = DummyEngine(save_commands=True)
        test_eng = MainEngine(backend=Simulator(),
                              engine_list=[AutoReplacer(rule_set)])

        correct_qb = correct_eng.allocate_qureg(n)
        correct_eng.flush()
        test_qb = test_eng.allocate_qureg(n)
        test_eng.flush()

        correct_eng.backend.set_wavefunction(basis_state, correct_qb)
        test_eng.backend.set_wavefunction(basis_state, test_qb)

        DenseHouseholderDec(gate_matrix) | test_qb
        test_gate | correct_qb

        test_eng.flush()
        correct_eng.flush()

        #assert correct_dummy_eng.received_commands[2].gate == test_gate
        #assert test_dummy_eng.received_commands[2].gate != test_gate

        for fstate in ['00', '01','10','11']:
            test = test_eng.backend.get_amplitude(fstate, test_qb)
            correct = correct_eng.backend.get_amplitude(fstate, correct_qb)
            assert correct == pytest.approx(test, rel=1e-12, abs=1e-12)

        All(Measure) | test_qb
        All(Measure) | correct_qb


@pytest.mark.parametrize("gate_matrix", create_iso(2,3))
def test_sparse_decomposition(gate_matrix):
    for basis_state in ([1, 0,0,0], [0,0,0, 1],[0,1,0,0],[0,0,1,0]):
        # Create single qubit gate with gate_matrix
        n=2
        test_gate = MatrixGate()
        test_gate.matrix = np.matrix(gate_matrix)

        correct_dummy_eng = DummyEngine(save_commands=True)
        correct_eng = MainEngine(backend=Simulator(),
                                 engine_list=[correct_dummy_eng])
        rule_set = DecompositionRuleSet(modules=[projectq.libs.math,
                                                 projectq.setups.decompositions])
        test_dummy_eng = DummyEngine(save_commands=True)
        test_eng = MainEngine(backend=Simulator(),
                              engine_list=[AutoReplacer(rule_set)])

        correct_qb = correct_eng.allocate_qureg(n)
        correct_eng.flush()
        test_qb = test_eng.allocate_qureg(n)
        test_eng.flush()

        correct_eng.backend.set_wavefunction(basis_state, correct_qb)
        test_eng.backend.set_wavefunction(basis_state, test_qb)

        SparseHouseholderDec(gate_matrix) | test_qb
        test_gate | correct_qb

        test_eng.flush()
        correct_eng.flush()

        #assert correct_dummy_eng.received_commands[2].gate == test_gate
        #assert test_dummy_eng.received_commands[2].gate != test_gate

        for fstate in ['00', '01','10','11']:
            test = test_eng.backend.get_amplitude(fstate, test_qb)
            correct = correct_eng.backend.get_amplitude(fstate, correct_qb)
            assert correct == pytest.approx(test, rel=1e-12, abs=1e-12)

        All(Measure) | test_qb
        All(Measure) | correct_qb