import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

n_ATOMS = 4000
filename = ""
###########################
#data = pd.read_table(filename, on_bad_lines="skip", comment="I", skiprows=9, delimiter=" ",
#                         names=["atom_nr", "type", "x", "y", "z", "vx", "vy", "vz", "fx", "fy", "fz", "pe"]).dropna()
q_Al = 1.4175
a_Al = 0.7852
b_Al = 0.034
c_Al = 36.82

q_O = -0.9450
a_O = 1.8215
b_O = 0.138
c_O = 90.61

coeff_Al = [q_Al, a_Al, b_Al, c_Al]
coeff_O = [q_O, a_O, b_O, c_O]

def matsui_pot(r, coeff_1, coeff_2):
    _epsilon = 4172.8356822474225 #kJ/(mol Å) 1/epsilon
    f = 4.184 #kJ/(mol Å)
    coloumb = coeff_1[0]*coeff_2[0]/(4*np.pi*r)*_epsilon

    vdw = coeff_1[3]*coeff_2[3]/r**6

    repulsion = f*(coeff_1[2] + coeff_2[2])*np.exp((coeff_1[1] + coeff_2[1] - r)/(coeff_1[2] + coeff_2[2]))
    #plt.figure()
    #plt.plot(r, coloumb, label="coloumb")
    #plt.plot(r, -vdw, label="vdw")
    #plt.plot(r, repulsion, label="repulsion")
    #plt.legend()
    #plt.show()

    return coloumb - vdw + repulsion


def plot_matsui(coeff_Al, coeff_O):
    r = np.arange(0.01, 7, 0.0001)
    pot_AlAl = matsui_pot(r, coeff_Al, coeff_Al)
    pot_AlO = matsui_pot(r, coeff_Al, coeff_O)
    pot_OO = matsui_pot(r, coeff_O, coeff_O)
    x_line = np.full(2, r[pot_AlAl == max(pot_AlAl)])
    y_line = [-0.4e9, 0.3e9]


    plt.figure()
    plt.title("Matsui potential", size=16)
    plt.plot(x_line, y_line, color="k", linestyle=":", linewidth=1.5)
    plt.plot(r, pot_AlAl, linewidth=2, color="teal", label="AlAl", alpha=0.8)
    plt.plot(r, pot_AlO, linewidth=2, color="darkblue", linestyle="--", label="AlO")
    plt.plot(r, pot_OO, linewidth=2, color="firebrick", linestyle="-.", label="OO")
    plt.xlabel("Interatomic distance [Å]")
    plt.ylabel("Energy [kJ/mol]")
    plt.xlim(0.8, 6)
    plt.ylim(-400, 600)
    #plt.ylim(-0.2e9, 0.3e9)
    plt.legend(frameon=False)
    plt.show()

    fig, (ax0, ax1, ax2) = plt.subplots(nrows=1, ncols=3, sharex=True, figsize=(12,6))

    ax0.set_title("AlAl")
    ax0.plot(r, 15/r**24, "--", label="short-range repulsive term")
    ax0.plot(r, pot_AlAl, label="Matsui potential")
    ax0.set_xlim(0.1, 3.5)
    ax0.set_ylim(-400, 600)
    ax0.legend()


    ax1.set_title("AlO")
    ax1.plot(r, 15/r**24, "--")
    ax1.plot(r, pot_AlO)
    ax1.set_xlim(0.1, 4.0)
    ax1.set_ylim(-400, 600)

    ax2.set_title("OO")
    ax2.plot(r, 15/r**24, "--")
    ax2.plot(r, pot_OO)
    ax2.set_xlim(0.1, 4.0)
    ax2.set_ylim(-400, 600)

    plt.show()


plot_matsui(coeff_Al, coeff_O)

