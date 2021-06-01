import sys
import numpy as np


"""
This script can be used test wheter the dump from two runs if exact paramters matches
Use case: testing code
"""


def main():
    if(len(sys.argv) < 2):
        print("file name missing\n")
        exit()
    # fname1 = '../dump/dump_08042020/fig7_1/dump_{}.txt'
    # fname2 = '../dump/dump_09042020/test_new_method/dump_{}.txt'
    # fname1 = fname1.format(int(sys.argv[1]))
    # fname2 = fname2.format(int(sys.argv[1]))
    fname1 = sys.argv[1]
    fname2 = sys.argv[2]
    d1=np.loadtxt(fname1)
    d2=np.loadtxt(fname2)
    truth=d1==d2
    idx = np.argwhere(truth==False)
    print("idx: ", idx)


if __name__ == '__main__':
    main()
