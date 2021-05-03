import projectq.libs.math
import projectq.setups.decompositions
from projectq.backends import Simulator, ResourceCounter
from projectq.cengines import (AutoReplacer, DecompositionRuleSet,
                               InstructionFilter, LocalOptimizer,
                               MainEngine, TagRemover)
from projectq.libs.math import (AddConstant, AddConstantModN,
                                MultiplyByConstantModN)
from projectq.meta import Control
from projectq.ops import (All, BasicMathGate, get_inverse, H, Measure, QFT, R,
                          Swap, X, DenseHouseholderDec, SparseHouseholderDec)
from projectq.meta._control import State
import time
from projectq.backends import Simulator, ResourceCounter
import numpy as np
    # make the compiler and run the circuit on the simulator backend


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

eng = MainEngine(Simulator(), compilerengines)
n = 2

x = eng.allocate_qureg(n)
iso = [[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]]
X | x[1]
DenseHouseholderDec(iso) | x




All(Measure) | x

eng.flush()
print('\n')
print('---------------------RESULTS--------------------------------')
for i in range(n):
    print('X_reg :',int(x[i]))

print(resource_counter)












#