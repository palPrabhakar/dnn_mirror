import sys
import numpy as np
import matplotlib.pyplot as plt
# from scipy.optimize import curve_fit
from scipy import optimize, special
# from 2d_peclet_number import get_reduced_d

def peclet_number_2d(x, omega):
    return np.sqrt(np.pi*x)*np.exp(x)*special.erfc(np.sqrt(x))-omega


def peclet_number(omega):
    sol = optimize.root_scalar(peclet_number_2d, args=omega, bracket=[0, 1], method='brentq')
    return sol.root

def get_reduced_d(omega):
    sol = peclet_number(omega)
    return 1/(2*sol)
    return a*np.power(x, b)

def main():
    if(len(sys.argv) < 3):
        print("Invalid argument.\nOption 1. -v <file_name>\nOption 2. -r <file_name>\n")
        exit(1)

    t, v, r, l = np.loadtxt(sys.argv[2], unpack=True)
    tp, vp = np.loadtxt('velocity.csv', unpack=True)
    tp, rp = np.loadtxt('radius.csv', unpack=True)
    t=t*1e-3

    if(sys.argv[1] == "-v"):
        plt.plot(t[1:], v[1:], 'r-', label='calculated')
        # plt.plot(tp, vp, 'b^', label='Tourret (2013)')
        # plt.title('Velocity, dx=10, omega=0.05, dt=1e-3')
        plt.ylabel('velocity')
    elif (sys.argv[1] == "-r"):
        plt.plot(t[1:], r[1:], 'r-', label='calculated')
        # plt.plot(tp, rp, 'b^', label='Tourret (2013)')
        # plt.title('Radius, dx=10, omega=0.05, dt=1e-3')
        plt.ylabel('radius')
    elif (sys.argv[1] == "-p"):
        print("Enter the omega value.\n")
        omega = float(sys.stdin.readline())
        # d=137.4469725 #omega=1
        # d=56.84926438 #omega=1.5
        # d=29.64150239 #omega=2
        d = get_reduced_d(omega)
        idx = 200
        vr = v*r;
        p = vr/(2*d)
        v_avg = np.average(v[idx:])
        r_avg = np.average(r[idx:])
        plt.plot(t[idx:], p[idx:], 'r-')
        plt.ylabel('Peclet number')
        print("Average Peclet number: ", np.average(p[idx:]))
        z = np.polyfit(t[idx:], p[idx:], 1)
        pfit = np.poly1d(z)
        plt.plot(t[idx:], pfit(t[idx:]))
        print("Average Peclet number: ", np.average(pfit(t[idx:])))
        print("Average values v, r: ", v_avg, r_avg)
        print("Average Peclet number: ", (v_avg*r_avg)/(2*d))
    else:
        # m=2000
        # n=-7000
        # dx=0.01
        # plt.plot(np.log10(t[m:n]), np.log10(l[m:n]*dx), '.')
        # pt, pc = curve_fit(fn, t[m:n], l[m:n]*dx)
        # print(pt)
        # plt.plot(np.log10(t[m:n]), np.log10(fn(t[m:n], *pt)))
        # plt.plot(t[m:n], l[m:n]*dx, '.')
        plt.plot(t[1:], l[1:], '.')
        plt.ylabel('len')


    plt.xlabel('time')
    # plt.legend()
    plt.show()


if __name__ == "__main__":
    main()



