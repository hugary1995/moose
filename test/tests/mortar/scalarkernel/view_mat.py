import pandas as pd
import numpy as np


def read_mat(filename):
    pbc = pd.read_csv(filename, index_col=False, sep="\s+")
    ii = pbc["i"].astype(int).to_numpy()
    jj = pbc["j"].astype(int).to_numpy()
    vv = pbc["v"].astype(np.float64).to_numpy()
    mat = np.zeros((7, 7))
    for idx, (i, j, v) in enumerate(zip(ii, jj, vv)):
        mat[i-1, j-1] = v
    return mat


np.set_printoptions(precision=3)
print("pbc")
print(read_mat("mat_pbc.csv"))
print("penalty pbc")
print(read_mat("mat_penalty_pbc.csv"))
