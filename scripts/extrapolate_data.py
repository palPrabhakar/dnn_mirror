import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

"""
To extrapolate the single needle radius, velocity values to long time step
Use case: Testing
"""

def func(x, a, b, c, d, e):
    return a*np.log(b*x) + c*np.exp(-d*x) + e


# def func(x, a, b, c):
#     return a*np.exp(-b*x) + c


# def func(x, a, b, c):
#     return a * np.log(b * x) + c


def main():
    if(len(sys.argv) < 2):
        print("Invalid Input\nFile name missing")
        exit(1)

    t, v, r, l = np.loadtxt(sys.argv[1], unpack=True)
    t=t*1e-4
    tt, vv = np.genfromtxt('paper_data_2.csv', delimiter=",", unpack=True)
    # print(len(tt))

    # popt, pcov = curve_fit(func, t[1:], v[1:])
    # tt = np.arange(0.1, 1000, 0.1)

    plt.plot(t[1:], v[1:], 'b.')
    # plt.plot(tt, func(tt, *popt), 'r-')
    plt.plot(tt[:-50], vv[:-50], 'r-')
    # plt.plot(t[1:], func(t[1:], *popt), 'r-')
    plt.grid()
    plt.axhline(1)
    plt.show()



if __name__ == '__main__':
    main()
