import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import convert_units


def data_to_file(x, y, dest):
    data = np.concatenate((x, y), axis=0)
    np.savetxt(fname=dest, X=data)

def read_data(fname):
    data = np.loadtxt(fname)
    rows = data.shape[0]
    x = data[:rows//2]
    y = data[rows//2:]

    return x, y



def plot_rdf_comparisons(f1=".\\rdf_data\\alal_avg_rdf_highT_Matsui.txt", f2=".\\rdf_data\\alal_avg_rdf_highT_Matsui.txt", f3=".\\rdf_data\\alal_avg_rdf_highT_DeepMD.txt", f4= ".\\rdf_data\\alal_avg_rdf_highT_DeepMD.txt"): #filenames

    r, lowT_Matsui = read_data(f1)
    _, highT_Matsui = read_data(f2)
    _, lowT_DeepMD = read_data(f3)
    _, highT_DeepMD = read_data(f4)

    plt.figure()
    plt.title("O-O rdf, T=2345K", size=16)
    plt.plot(r, highT_Matsui, "-o", color="r", linewidth=2, label="Matsui")
    plt.plot(r, highT_DeepMD, "-", color="k", linewidth=2, label="DeepMD")
    plt.xlabel("r [Å]", size=14)
    plt.ylabel("g(r)", size=14)
    plt.xlim(0, 9.5)
    plt.legend(fontsize=14)
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.show()


    plt.figure()
    plt.title("O-O rdf, T=300K", size=16)
    plt.plot(r, lowT_Matsui, "-o", color="r", linewidth=2, label="Matsui")
    plt.plot(r, lowT_DeepMD, "-", color="k", linewidth=2, label="DeepMD")
    plt.xlabel("r [Å]", size=14)
    plt.ylabel("g(r)", size=14)
    plt.xlim(0, 9.5)
    plt.legend(fontsize=14)
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.show()



def plot_adf_comparisons(f1=".\\adf_data\\alalal_avg_adf_lowT_Matsui.txt", f2=".\\adf_data\\alalal_avg_adf_highT_Matsui.txt", f3=".\\adf_data\\alalal_avg_adf_lowT_DeepMD.txt", f4= ".\\adf_data\\alalal_avg_adf_lowT_DeepMD.txt", title="O-O-O"): #filenames

    r, lowT_Matsui = read_data(f1)
    _, highT_Matsui = read_data(f2)
    _, lowT_DeepMD = read_data(f3)
    _, highT_DeepMD = read_data(f4)

    plt.figure()
    plt.title(title + " adf, T=2345K", size=16)
    plt.plot(r, highT_Matsui, "-o", color="r", linewidth=2, label="Matsui")
    plt.plot(r, highT_DeepMD, "-", color="k", linewidth=2, label="DeepMD")
    plt.xlabel("Degrees", size=14)
    plt.ylabel("Count", size=14)
    plt.legend(fontsize=14)
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.show()


    plt.figure()
    plt.title(title + " adf, T=300K", size=16)
    plt.plot(r, lowT_Matsui, "-o", color="r", linewidth=2, label="Matsui")
    plt.plot(r, lowT_DeepMD, "-", color="k", linewidth=2, label="DeepMD")
    plt.xlabel("Degrees", size=14)
    plt.ylabel("Count", size=14)
    plt.legend(fontsize=14)
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.show()


plot_adf_comparisons(f1=".\\adf_data\\alalal_avg_adf_lowT_Matsui.txt", f2=".\\adf_data\\alalal_avg_adf_highT_Matsui.txt", f3=".\\adf_data\\alalal_avg_adf_lowT_DeepMD.txt", f4= ".\\adf_data\\alalal_avg_adf_highT_DeepMD.txt", title="Al-Al-Al")
plot_adf_comparisons(f1=".\\adf_data\\aloal_avg_adf_lowT_Matsui.txt", f2=".\\adf_data\\aloal_avg_adf_highT_Matsui.txt", f3=".\\adf_data\\aloal_avg_adf_lowT_DeepMD.txt", f4= ".\\adf_data\\aloal_avg_adf_highT_DeepMD.txt", title="Al-O-Al")
plot_adf_comparisons(f1=".\\adf_data\\oalo_avg_adf_lowT_Matsui.txt", f2=".\\adf_data\\oalo_avg_adf_highT_Matsui.txt", f3=".\\adf_data\\oalo_avg_adf_lowT_DeepMD.txt", f4= ".\\adf_data\\oalo_avg_adf_highT_DeepMD.txt", title="O-Al-O")
plot_adf_comparisons(f1=".\\adf_data\\ooo_avg_adf_lowT_Matsui.txt", f2=".\\adf_data\\ooo_avg_adf_highT_Matsui.txt", f3=".\\adf_data\\ooo_avg_adf_lowT_DeepMD.txt", f4= ".\\adf_data\\ooo_avg_adf_highT_DeepMD.txt", title="O-O-O")


#plot_rdf_comparisons(f1=".\\rdf_data\\oo_avg_rdf_lowT_Matsui.txt", f2=".\\rdf_data\\oo_avg_rdf_highT_Matsui.txt", f3=".\\rdf_data\\oo_avg_rdf_lowT_DeepMD.txt", f4= ".\\rdf_data\\oo_avg_rdf_highT_DeepMD.txt")

def get_min(file, r_min, r_max):
    r, g = read_data(file)
    g_min = max(g[np.logical_and(r>r_min, r<r_max)])
    r_gmin = r[g == g_min]

    return r_gmin, g_min

def calc_coord_nr(file, rmax, n_atoms, l=31.82):
    r, g = read_data(file)
    rbz = r[r<rmax]
    g_peak = g[r<rmax]

    bz = np.mean(np.diff(r)) #approx 0.05Å

    coord_nr = np.sum(g_peak*rbz**2)

    return coord_nr*bz*4*np.pi*n_atoms/(l**3)


#rm, gm = get_min(file=".\\rdf_data\\oo_avg_rdf_lowT_Matsui.txt", r_min=2, r_max=3.5)
nAl = 1600
nO = 2400
l = 31.82

"""

print("\nAl-O")
print(calc_coord_nr(file=".\\rdf_data\\alo_avg_rdf_lowT_Matsui.txt", rmax=2.44, n_atoms=nO))
print(calc_coord_nr(file=".\\rdf_data\\alo_avg_rdf_lowT_Matsui.txt", rmax=2.44, n_atoms=nAl))

print(calc_coord_nr(file=".\\rdf_data\\alo_avg_rdf_highT_Matsui.txt", rmax=2.54, n_atoms=nO))
print(calc_coord_nr(file=".\\rdf_data\\alo_avg_rdf_lowT_DeepMD.txt", rmax=2.49, n_atoms=nO))
print(calc_coord_nr(file=".\\rdf_data\\alo_avg_rdf_highT_DeepMD.txt", rmax=2.64, n_atoms=nO))

print("\nAl-Al")
print(calc_coord_nr(file=".\\rdf_data\\alal_avg_rdf_lowT_Matsui.txt", rmax=3.93, n_atoms=nAl))
print(calc_coord_nr(file=".\\rdf_data\\alal_avg_rdf_highT_Matsui.txt", rmax=4.03, n_atoms=nAl))
print(calc_coord_nr(file=".\\rdf_data\\alal_avg_rdf_lowT_DeepMD.txt", rmax=3.88, n_atoms=nAl))
print(calc_coord_nr(file=".\\rdf_data\\alal_avg_rdf_highT_DeepMD.txt", rmax=4.08, n_atoms=nAl))

print("\nO-O")
print(calc_coord_nr(file=".\\rdf_data\\oo_avg_rdf_lowT_Matsui.txt", rmax=3.33, n_atoms=nO))
print(calc_coord_nr(file=".\\rdf_data\\oo_avg_rdf_highT_Matsui.txt", rmax=3.23, n_atoms=nO))
print(calc_coord_nr(file=".\\rdf_data\\oo_avg_rdf_lowT_DeepMD.txt", rmax=3.28, n_atoms=nO))
print(calc_coord_nr(file=".\\rdf_data\\oo_avg_rdf_highT_DeepMD.txt", rmax=3.28, n_atoms=nO))
"""