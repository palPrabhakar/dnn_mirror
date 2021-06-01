import sys
import matplotlib.pyplot as plt


def read_data(file_name):
    vel_arr = []
    rad_arr = []
    with open(file_name) as f:
        for line in f:
            d = line.split()
            if len(d) == 6:
                vel_arr.append(float(d[1].rstrip(',')))
                rad_arr.append(float(d[3].rstrip(',')))
    return vel_arr, rad_arr


def main():
    if (len(sys.argv) < 2):
        print("Input file name")
        exit(1)

    vel, rad = read_data(sys.argv[1])
    x = range(len(vel))
    # print(vel, rad)
    plt.plot(x[:-40], vel[:-40])
    # plt.plot(x[:-40], rad[:-40])
    plt.show()


if __name__ == '__main__':
    main()
