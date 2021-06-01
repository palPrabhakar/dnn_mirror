import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize, special

def peclet_number_2d(x, omega):
    return np.sqrt(np.pi*x)*np.exp(x)*special.erfc(np.sqrt(x))-omega


def peclet_number(omega):
    sol = optimize.root_scalar(peclet_number_2d, args=omega, bracket=[0, 1], method='brentq')
    return sol.root

def main():
    omega = [0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    peclet = [0.003587035, 0.0083690187, 0.015811866, .040591157, .082328875, 0.131043339, 0.226224843, .4127633395]
    peclet_th = []
    rel_err = []

    for val in omega:
        peclet_th.append(peclet_number(val))

    for i in range(len(omega)):
        err = (np.abs(peclet_th[i] - peclet[i])/peclet_th[i])*100
        rel_err.append(err)

    print(rel_err)
    plt.plot(omega, rel_err, 'o', fillstyle='none')
    plt.xlabel('Supersaturation')
    plt.ylabel('Relative Error (%)')
    plt.show()


if __name__ == '__main__':
    main()

