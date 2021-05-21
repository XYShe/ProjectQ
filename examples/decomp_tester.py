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
from projectq.meta._control import State
import time
from projectq.backends import Simulator, ResourceCounter, CircuitDrawerMatplotlib
import numpy as np
from projectq.libs._utils import chop,ctrlmat,stdHHref
    # make the compiler and run the circuit on the simulator backend
from scipy.stats import unitary_group

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
x = eng.allocate_qureg(n)
iso = [[Complex(-0.1742226972380019,0.22676721299188038),
        Complex(0.30431307206823016,-0.7345582542754363),
        Complex(0.18555117703730004,0.13799045645485278),
        Complex(0.18490101079089102,-0.4454007397380378)],
       [Complex(0.181062856534321,-0.0649957253285508),
        Complex(0.40851517573037993,0.036355752826448866),
        Complex(0.28033960283698156,-0.8253581014523761),
        Complex(-0.16316366736793964,-0.09141519133176651)],
       [Complex(-0.6019898770916067,0.0827892509304244),
        Complex(0.21719916999719355,0.08459761540589042),
        Complex(0.04209922828446449,0.1305498526254742),
        Complex(-0.7459175716760235,0.034831868339822186)],
       [Complex(-0.6935457294596152,-0.1759640006069421),
        Complex(-0.37608263965928074,-0.061975380964489436),
        Complex(0.21637369278080254,-0.34797977255194745),
        Complex(0.38281878122410984,0.16820872643571408)]]
#X | x[0]
with Control(eng,x[0]):
    X | x[1]
#SparseHouseholderDec(iso) | x

#StatePreparation([0,1/np.sqrt(2),1/np.sqrt(2),0]) | x

All(Measure) | x
print(eng.backend.allgate())

eng.flush()
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


#