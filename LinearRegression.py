import pandas as pd
import numpy as np
from statsmodels.regression.linear_model import OLS

# H_set = np.loadtxt("H1.csv", delimiter=",")
Z_set = np.loadtxt("Z1.csv", delimiter=",")
A_set = np.loadtxt("A1.csv", delimiter=",")
Y_set = np.loadtxt("Y1.csv", delimiter=",")

ZA_stack = np.vstack([Z_set, A_set, np.ones(Z_set.shape[1])])

reg = OLS(Y_set, np.transpose(ZA_stack), hasconst = True).fit()
ATE = reg.params[-2]
print(ATE)
