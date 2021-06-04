import projectq.libs.math
import projectq.setups.decompositions
from projectq.backends import Simulator, ResourceCounter
from projectq.cengines import (AutoReplacer, DecompositionRuleSet,
                               InstructionFilter, LocalOptimizer,
                               MainEngine, TagRemover)
from projectq.libs.math import (AddConstant, AddConstantModN,
                                MultiplyByConstantModN)
from projectq.meta import Control
from projectq.ops import (All, BasicMathGate, get_inverse, H, Measure, QFT, R, DaggeredGate,
                          Swap, X, DenseHouseholderDec, SparseHouseholderDec, StatePreparation,C,
                          CNOT,MatrixGate)
import time
import random
from projectq.backends import Simulator, ResourceCounter, CircuitDrawerMatplotlib
import numpy as np
from projectq.libs._utils import chop,ctrlmat,stdHHref, create_iso
    # make the compiler and run the circuit on the simulator backend
from scipy.stats import unitary_group
import scipy
resource_counter = ResourceCounter()
rule_set = DecompositionRuleSet(modules=[projectq.libs.math,
                                         projectq.setups.decompositions])


compilerengines = [AutoReplacer(rule_set),
                   TagRemover(),
                   LocalOptimizer(3),
                   AutoReplacer(rule_set),
                   TagRemover(),
                   LocalOptimizer(3),
                   resource_counter]
def Complex(a,b):
    return a+1j*b
eng = MainEngine(Simulator())
n = 2
x = eng.allocate_qureg(1)
y = eng.allocate_qureg(1)

iso = [[0.395249 - 0.527105*1j, 0, 0.746712 - 0.0914316*1j,0],
[0, -0.106333 + 0.994331*1j, 0, 0],
[-0.0130445 + 0.0312883 *1j, 0,0.0276758 - 0.0107422 *1j, -0.962676 + 0.26688 *1j],
[-0.480347 - 0.577976 *1j, 0, -0.0736539 + 0.65403 *1j, -0.00947929 + 0.0440523 *1j]]
#X | x[0]
iso1 = [[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]]
x.extend(y)
print(x)
DenseHouseholderDec(iso1) | x

#StatePreparation([0,1/np.sqrt(2),1/np.sqrt(2),0]) | x


eng.flush()
mat=  eng.backend.matout()
chop(mat)
print(mat)

All(Measure) | x
#st = np.array(eng.backend.cheat()[-1])
#chop(st)

print('\n')
print('---------------------RESULTS--------------------------------')
#for i in range(2**n):
#    print(st[i])

for i in range(n):
    print(int(x[i]))

#print(resource_counter)

#print(eng.backend.allgate())

#print(ctrlmat([1],[0,2],3,np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])))

#vec = np.array([[-0.48852773+0.63586475j,  0.11290145-0.04052798j, -0.37536981+0.0516231j,
     #  -0.43245931-0.10972207j]])

#print(stdHHref(vec.transpose(),np.array(iso)))


iso = create_iso(2,3)
print(iso)