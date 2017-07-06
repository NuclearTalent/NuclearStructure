import numpy as np

def binomial(n, k):
    uvector = np.array(range(k+1, n+1))
    lvector = np.array(range(1, n-k+1))

    return np.e**(np.log(uvector).sum() - np.log(lvector).sum())

print binomial(200, 20)**2
