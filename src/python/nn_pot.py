import numpy as np
import matplotlib.pyplot as plt
#import structure_analysis
#import lammps_utils

###########################

def weight_func(rc, rcs):
    r1 = np.arange(0, rcs, 0.1)
    r2 = np.arange(rcs, rc, 0.1)
    return r1, 1/r1, r2, 1/r2*(0.5*np.cos(np.pi*(r2 - rcs)/(rc - rcs)) + 0.5)

r_c = [6]
r_cs = [0.5, 1, 1.5, 2, 2.5, 3, 3.5]

plt.figure()
for rc in r_c:
    for rcs in r_cs:
        r1, s1, r2, s2 = weight_func(rc, rcs)
        r = np.concatenate((r1, r2))
        s = np.concatenate((s1, s2))
        plt.plot(r, s, label ="$r_{cs}$ = "+ str(rcs))
plt.legend(fontsize=16)
plt.title("Weight function with different $r_{cs}$", size=18)
plt.ylabel("Weight", size=16)
plt.xlabel("r [Å]", size=16)
plt.xlim(1, 5.5)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim(0, 1)
plt.show()

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

    return coloumb - vdw + repulsion


def plot_matsui(coeff_Al, coeff_O):
    r = np.arange(0.01, 7, 0.0001)
    pot_AlAl = matsui_pot(r, coeff_Al, coeff_Al)
    pot_AlO = matsui_pot(r, coeff_Al, coeff_O)
    pot_OO = matsui_pot(r, coeff_O, coeff_O)
    x_line = np.full(2, r[pot_AlAl == max(pot_AlAl)])
    y_line = [-0.4e9, 0.3e9]


    plt.figure(figsize=(6,4))
    plt.title("Matsui potential", size=18)
    #plt.plot(x_line, y_line, color="k", linestyle=":", linewidth=1.5)
    plt.plot(r, pot_AlAl, linewidth=2, color="teal", label="Al-Al", alpha=0.8)
    plt.plot(r, pot_AlO, linewidth=2, color="darkblue", linestyle="--", label="Al-O")
    plt.plot(r, pot_OO, linewidth=2, color="firebrick", linestyle="-.", label="O-O")
    plt.xlabel("Interatomic distance [Å]", size=16)
    plt.ylabel("Energy [kJ/mol]", size=16)
    plt.xlim(0.01, 3)
    #plt.ylim(-400, 600)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylim(-5e3, 2e4)
    plt.legend(frameon=False, fontsize=16)
    plt.show()

    fig, (ax0, ax1, ax2) = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, figsize=(12,6))
    ax0.tick_params(axis="both", labelsize=14)
    ax0.set_ylabel("Energy [kJ/mol]", size=16)
    ax0.set_title("Al-Al", size=18)
    ax0.plot(r, 15/r**24, color="red", linestyle="--", label="Repulsive term")
    ax0.plot(r, pot_AlAl, color="k", label="Matsui potential")
    ax0.set_xlim(0.1, 3.5)
    ax0.set_ylim(-400, 600)
    ax0.legend(fontsize=16)



    ax1.set_title("Al-O", size=18)
    ax1.plot(r, 15/r**24, color="red", linestyle="--", label="Repulsive term")
    ax1.plot(r, pot_AlO, "k")
    ax1.tick_params(axis="both", labelsize=14)
    ax1.set_xlim(0.1, 4.0)
    ax1.set_ylim(-400, 600)
    print("Minima: ", r[r>1][pot_AlO[r>1]==min(pot_AlO[r>1])])

    ax2.set_title("O-O", size=18)
    ax2.tick_params(axis="both", labelsize=14)
    ax2.plot(r, 15/r**24, color="red", linestyle="--", label="Repulsive term")
    ax2.plot(r, pot_OO, "k")
    ax2.set_xlim(0.1, 4.0)
    ax2.set_ylim(-400, 600)

    plt.show()


plot_matsui(coeff_Al, coeff_O)




fname = "C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\data.dump"
r, types_n, f, e = structure_analysis.get_dump(fname, n_atoms=2, return_all=True)

n_frames = r.shape[1]

dist = np.sqrt(np.sum((r[0, :] - r[1, :])**2, axis=-1))
CU = lammps_utils.ConvertUnits()
print(dist[e[0]+e[1] == np.min(e[0, :]+e[1, :])])
m_pot = CU.kcal_mol_to_eV(CU.kJ_to_kcal(matsui_pot(dist, coeff_Al, coeff_O)))
print(min(m_pot[dist>1]))
plt.figure()
plt.title("NN pot and Matsui plot, Al-O")
plt.plot(dist, (e[0, :]+e[1, :])/2 + (min(m_pot[dist>1]) - min((e[0]+e[1])[dist>1])/2), "-o", label="Pot_E avg")
#plt.plot(dist, e[0], "-o", label="Pot_E, Al")
#plt.plot(dist, e[1], "-o", label="Pot_E, O")
plt.plot(dist, m_pot, label="Matsui")
#plt.plot(np.full(len(dist), dist[e[0]+e[1] == np.min(e[0, :]+e[1, :])]), np.linspace(-30, 7, num=len(dist))) #line through minima
plt.ylim(-5, 5)
plt.xlim(0.8, 8)
plt.ylabel("Energy [eV]")
plt.xlabel("Distance [Å]")
plt.grid()
plt.legend()
plt.show()

def get_force(d, e):
    force = np.zeros(len(d)-1)
    print(d.shape)
    for i in range(len(force)-1):
        if d[i] - d[i-1] > 0:
            force[i] = -(e[i+1] - e[i])/(d[i+1]-d[i])
        else:
            force[i] = force[i-1]
    return force

dist_m = np.arange(0.5, 8, 0.1)
print(dist_m.shape)
matsui_p = CU.kcal_mol_to_eV(CU.kJ_to_kcal(matsui_pot(dist_m, coeff_Al, coeff_O)))

m_force = get_force(dist_m, matsui_p)

plt.figure()
plt.title("Force between Al and O")
#plt.plot(dist, f[0, :, -1], "-o", label="Al")
plt.plot(dist, 2*f[1, :, -1], "-o", label="O")
plt.plot(dist_m[:-1], m_force, label="Matsui")
#plt.plot(dist, CU.kcal_mol_to_eV(CU.kJ_to_kcal(matsui_pot(dist, coeff_Al, coeff_O))), label="Matsui")
#plt.plot(np.full(len(dist), dist[e[0]+e[1] == np.min(e[0, :]+e[1, :])]), np.linspace(-30, 7, num=len(dist))) #line through minima
plt.ylim(-5, 90)
plt.xlim(0.5, 3.5)
plt.ylabel("Force [eV/Å]")
plt.xlabel("Distance [Å]")
plt.grid()
plt.legend()
plt.show()

def get_pot(dist, f):
    pot = np.zeros(len(dist)-1)
    for i in range(1, len(pot)):
        pot[i] = pot[i-1] - 2*f[0, i-1, -1]*np.abs((dist[i]-dist[i-1]))

    return pot



pot = get_pot(dist, f)

plt.figure()
plt.title("Potential between Al and O")
plt.plot(dist, m_pot, label="Matsui")
plt.plot(dist[:-1], pot, label="NN Pot")
#plt.plot(dist[:-1], -f[0, :-1, -1]*np.diff(dist), "-o", label="Al-O")
#plt.plot(dist, f[1, :, -1]*dist, "-o", label="O")
#plt.plot(dist, CU.kcal_mol_to_eV(CU.kJ_to_kcal(matsui_pot(dist, coeff_Al, coeff_O))), label="Matsui")
#plt.plot(np.full(len(dist), dist[e[0]+e[1] == np.min(e[0, :]+e[1, :])]), np.linspace(-30, 7, num=len(dist))) #line through minima
plt.ylim(-5, 5)
plt.xlim(1, 8)
plt.ylabel("Energy [eV]")
plt.xlabel("Distance [Å]")
plt.grid()
plt.legend()
plt.show()

