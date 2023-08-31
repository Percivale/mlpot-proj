import numpy as np
import matplotlib.pyplot as plt
import structure_analysis as sa

def data_to_file(x, y, dest):
    data = np.concatenate((x, y), axis=0)
    np.savetxt(fname=dest, X=data)

def read_data(fname):
    data = np.loadtxt(fname)
    x = data[:, 0]
    y = data[:, 1]

    return x, y

def plot_rdf_comparisons(f1, f2, f3, f4): #filenames

    lowT_Matsui, r = read_data(f1)
    highT_Matsui, r = read_data(f2)
    lowT_DeepMD, r = read_data(f3)
    highT_DeepMD, r = read_data(f4)

    plt.figure()
    plt.title("Avg rdf over 20 ps, T=2345K")
    plt.plot(r, highT_Matsui, label="Matsui")
    plt.plot(r, highT_DeepMD, label="DeepMD")
    plt.show()

    plt.figure()
    plt.title("Avg rdf over 20ps, T=300K")
    plt.plot(r, lowT_Matsui, label="Matsui")
    plt.plot(r, lowT_DeepMD, label="DeepMD")
    plt.show()

    return 0