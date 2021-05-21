import numpy as np
from projectq.ops import Toffoli, Rx, Rz, Ry, Reflection
import itertools


def normalized(a, axis=-1, order=2):
    l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
    l2[l2==0] = 1
    return a / np.expand_dims(l2, axis)

def schmidtDec(v):

    v = np.array(v)


    state_dims = v.shape
    v = v.reshape((state_dims[0]//2,state_dims[0]//2))

    mindim = np.min(v.shape)

    vecs1, diags, vecs2_h = np.linalg.svd(v)
    vecs2 = vecs2_h.transpose()
    decomposition = [(diags[k], vecs1[:, k], vecs2[:, k])
                     for k in range(mindim)]

    decomposition = sorted(decomposition, key=lambda dec: dec[0], reverse=True)

    return decomposition

def ComputeHouseholderVec(v,i):
    v2 = v
    v2[i] -= np.exp(-1j*(np.pi - np.angle(v[i]) ))

    v2 = v2/np.linalg.norm(v2,ord=2)
    return v2


def nextGen(V):
    m = V.shape[0]
    n = V.shape[1]
    cols = np.array([V[:,i].getnnz() for i in range(n)])


    rows  = np.array([V[i,:].getnnz() for i in range(m)])

    if np.any(cols>1):
        mincol = np.amin(cols[cols>1])

        j = np.where(cols==mincol)[0][0]
    else:
        j = 0

    if np.any(rows>1):

        minrow = np.amin(rows[rows>1])
        i = np.where(rows==minrow)[0][0]
    else:
        i = 0
    print(i,j)
    return (i,j)


def chop(v, tol=10**(-8)):
    v.real[abs(v.real) < tol] = 0.0
    v.imag[abs(v.imag) < tol] = 0.0

def stdHHref(vec,V):

    result = V - 2*(np.dot(vec,np.dot(vec.conj().T,V)))
    chop(result)
    return result

def isPermutedDiag(V):
    coln = V.shape[-1]
    nzeros = [V[:,i].getnnz() for i in range(coln)]
    return max(nzeros) == 1

def applyGatetoISO(gates,iso):

    for mat in gates:
        iso = np.matmul(mat,iso)
    return iso

def ctrlmat(ctrl,target,n, gate = np.array([[0,1],[1,0]])):
    gatedict = {}
    gatedict[0]= np.array([[1,0],[0,0]])
    gatedict[1] = np.array([[0,0],[0,1]])
    gatedict[-1] = np.array([[1,0],[0,1]])
    gatedict[-2] = gate

    nctrl = len(ctrl)
    perms = np.ones((2**nctrl,n)) * -1

    lst = list(itertools.product([0, 1], repeat=nctrl))
    perms[:,ctrl] = lst
    perms[-1,target] = -2
    prod = 0
    for i in range(perms.shape[0]):
        row = gatedict[perms[i,0]]
        for j in range(1,perms.shape[1]):
            row = np.kron(row,gatedict[perms[i,j]])
        prod += row
    return prod

def decBitonic(i,j,list):
    d=1<<(i-j)
    n=int(np.log2(len(list)))
    gates = []
    for k in range(0,2**n):
        up = ((k>>i) & 2) == 0
        a = k+1
        b = (k | d) + 1
        if (k & d) == 0:
            if (list[a-1] > list[b-1]) == up:
                list[a-1], list[b-1]= list[b-1], list[a-1]
                target = n-i+j-1
                ctrl_qubits = np.setdiff1d(np.array(range(n)),[target])
                prod = ctrlmat(ctrl_qubits,target,n)
                gates.append(prod)
    return gates

def decPermDiag(V):
    n=int(np.log2(V.shape[0]))
    m=int(np.log2(V.shape[1]))
    gates = []
    nzeropos = [V[:,i].nonzero()[0][0] for i in range(2**m)]
    if nzeropos == list(np.sort(nzeropos)):
        return []
    else:
        for i in range(m):
            for j in range(i+1):
                gates.extend( decBitonic(i,j,nzeropos))

    return gates
