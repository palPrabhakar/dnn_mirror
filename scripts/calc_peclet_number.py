import sys
import numpy as np
import matplotlib.pyplot as plt

"""
This script is absolete not now
the fucntionality merged in plot_vel_rad.py
"""

def main():
    if len(sys.argv) < 2:
        print("Wrong, input file name missing.\n")
        exit(1)

    t, v, r, l = np.loadtxt(sys.argv[1], unpack=True)
    t=t*1e-2
    omega=0.05
    d=np.pi/(2*omega**2)
    print(d)
    print(1/(2*d))
    print(len(t))
    peclet_number=np.array(v*r/(2*d))
    print("peclet_number: ", np.average(peclet_number[50:]))
    plt.plot(t, peclet_number)
    plt.show()


if __name__ == "__main__":
    main()

