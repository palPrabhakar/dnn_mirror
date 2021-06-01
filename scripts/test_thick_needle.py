import numpy as np
import sys

"""
This script is used to check if the thickness of the needle
in the dump output matches that calculated using the growth conditoin.
This ensures that right needle BC is used in the transient solver
"""

def weak_test(file_name, length, radius, yidx):
    d = np.loadtxt(file_name)
    for i in range(length+1):
        wd = int(np.sqrt(2*radius*(length-i)))
        if d[yidx+wd+1, i] == 0:
            return False
    return True


def strong_test(file_name, length, radius, yidx):
    # print(file_name, length, radius)
    d = np.loadtxt(file_name)
    needle_index_calc = []
    for i in range(length+1):
        wd = int(np.sqrt(2*radius*(length-i)))
        needle_index_calc.append([yidx, i])
        for j in range(1, wd+1):
            needle_index_calc.append([yidx+j, i])
            needle_index_calc.append([yidx-j, i])
    zero_index = np.argwhere(d==0)
    # print("calc: ", needle_index_calc.sort())
    # print("zero_index: ", zero_index)
    needle_index_calc.sort()
    # truth1 = zero_index == needle_index_calc
    # print("zero_index: ", len(zero_index), "calc: ", len(needle_index_calc))
    truth3 = len(zero_index) == len(needle_index_calc)
    truth2 = np.array_equal(zero_index, needle_index_calc)

    return truth2 and truth3


def main():
    if(len(sys.argv) < 3):
        print("Wrong input\nEnter file name and y position of needle.\n")

    str_dir = sys.argv[1]
    # print(str_dir)
    data_file = str_dir + '/dump_vel_rad.dat'
    _iter, vel, rad, length = np.loadtxt(data_file, unpack=True)

    for i in range(len(_iter)):
        file_name = str_dir + '/dump_' + str(int(_iter[i])) + '.txt'
        # truth = weak_test(file_name, int(length[i]), rad[i], int(sys.argv[2]))
        truth = strong_test(file_name, int(length[i]), rad[i], int(sys.argv[2]))
        print(_iter[i], truth)


if __name__ == '__main__':
    main()


