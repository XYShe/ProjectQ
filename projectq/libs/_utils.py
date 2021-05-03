import numpy as np

def normalized(a, axis=-1, order=2):
    l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
    l2[l2==0] = 1
    return a / np.expand_dims(l2, axis)

def SchmidtDec(v):
    v = np.array(v)


    state_dims = v.shape
    v = v.reshape((state_dims[0]//2,state_dims[0]//2))
    mindim = np.min(v.shape)

    vecs1, diags, vecs2_h = np.linalg.svd(v)
    vecs2 = vecs2_h.transpose()
    print(vecs1, diags, vecs2)
    decomposition = [(diags[k], vecs1[:, k], vecs2[:, k])
                     for k in range(mindim)]

    decomposition = sorted(decomposition, key=lambda dec: dec[0], reverse=True)

    return decomposition

def ComputeHouseholderVec(v,i):
    v2 = v
    v2[i] -= np.exp(-1j*(np.pi - np.angle(v[i]) ))

    v2 = v2/np.linalg.norm(v2,ord=2)
    return v2

