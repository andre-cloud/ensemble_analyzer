#%%
import matplotlib.pyplot as plt
import numpy as np

X = np.linspace(1, 9, 10**4)
LIM_EXP1, LIM_EXP2 = 1.3, 8.5
INT1, INT2 = 2, 5.5

Y = np.zeros_like(X)

def gau(x0, sigma, X):
    return np.exp(-0.5*((X - x0)**2)/sigma**2)

def sigma_from_value(x1, x0, v, A=1.0):
    d = abs(x1 - x0)
    return d / np.sqrt(2 * np.log(A/v))

# sigma per code che valgono 0.2 ai limiti
sigma1 = sigma_from_value(LIM_EXP1, INT1, v=0.2, A=1)
sigma2 = sigma_from_value(LIM_EXP2, INT2, v=0.2, A=1)

gau1 = gau(INT1, sigma1, X)
gau2 = gau(INT2, sigma2, X)

Y[(LIM_EXP1 < X) & (X < INT1)] = gau1[(LIM_EXP1 < X) & (X < INT1)]
Y[(INT1 <= X) & (X <= INT2)] = 1
Y[(INT2 < X) & (X < LIM_EXP2)] = gau2[(INT2 < X) & (X < LIM_EXP2)]

Y[X < LIM_EXP1] = 0
Y[X > LIM_EXP2] = 0

plt.plot(X, Y)
plt.show()
# %%
