import sys
import matplotlib.pyplot as plt
import numpy as np

"""
To plot the concentration contour
"""


def plot_contour(file_name):
    fig, ax = plt.subplots()
    file_name = sys.argv[1]
    data = np.loadtxt(file_name)
    print("max value: ", np.amax(data))
    # print(data[200, data.shape[1]-1])
    x = np.arange(0, data.shape[1], 1)
    y = np.arange(0, data.shape[0], 1)
    X, Y = np.meshgrid(x, y)
    # cs = plt.contour(X, Y, data, np.arange(0, 0.006, 0.001))
    # cs = plt.pcolormesh(X, Y, data)
    # cs = plt.contour(X, Y, data, [0, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0], colors=('green', 'blue', 'magenta', 'cyan', 'red', 'yellow', 'purple', 'orange', 'brown', 'black'))
    # cs = plt.contour(X, Y, data, [0, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0])
    # ax.clabel(cs, inline=1, fontsize=10)
    plt.contourf(data, 50)
    # plt.imshow(data)
    plt.gca().set_aspect('equal')
    plt.set_cmap('jet')
    plt.colorbar()
    # plt.savefig('contour_plot.png')
    plt.grid()
    plt.show()


# def plot_FIF():
#     file_name = './dump_ss_1m_iterm/flux_{}.dat'
#     file_names = ['h3', 'h5', 'h7', 'h10']
#     # file_names = ['h3']
#     rel_err = []
#     fth = np.sqrt((30*0.3236*0.3236)/np.pi)
#     print("fth", fth)
#     for name in file_names:
#         fname = file_name.format(name)
#         num, val = np.loadtxt(fname, unpack=True)
#         favg = np.average(val[-10:])
#         fcal = np.sqrt(np.average(val[-10:]))
#         rel_err.append(np.abs(fth - fcal)/fth)
#         # print(name, rel_err[name], fcal, fname, favg)
#     x = [24, 40, 56, 80]

#     plt.plot(x, rel_err, 'x')
#     plt.show()
#     # print(data['h3'])


def main():
    if len(sys.argv) < 2:
        print("Input valid file name")
        exit(1)

    plot_contour(sys.argv[1])
    # plot_FIF();


if __name__ == '__main__':
    main()
