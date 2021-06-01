import numpy as np
from scipy import optimize, special
import sys


def peclet_number_2d(x, omega):
    return np.sqrt(np.pi*x)*np.exp(x)*special.erfc(np.sqrt(x))-omega


def peclet_number(omega):
    sol = optimize.root_scalar(peclet_number_2d, args=omega, bracket=[0, 1], method='brentq')
    return sol.root

def get_reduced_d(omega):
    sol = peclet_number(omega)
    return 1/(2*sol)


def main():
    if(len(sys.argv) < 2):
        print("Enter supersaturation value.\n")
        exit()

    ss = float(sys.argv[1])
    print(peclet_number(ss))
    print(get_reduced_d(ss))


if __name__ == '__main__':
    main()


