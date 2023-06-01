import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def get_data(filename: str):
    data = pd.read_table(filename, on_bad_lines="skip", comment="I", skiprows=9, delimiter=" ",
                         names=["atom_nr", "type", "x", "y", "z"]).dropna()  # .reset_index(drop=True)

    x = np.asarray(data["x"])
    y = np.asarray(data["y"])
    z = np.asarray(data["z"])
    types = np.asarray(data["type"], dtype=int)
    atom_idx = np.asarray(data["atom_nr"], dtype=int)

    n_frames = int(len(x) // 4000)
    r = np.zeros((4000, n_frames, 3))
    types_n = np.zeros((4000, n_frames), dtype=int)
    atom_idx_n = np.zeros((4000, n_frames), dtype=int)

    for i in range(n_frames):
        r[:, i, 0] = x[i * 4000:4000 + i * 4000]
        r[:, i, 1] = y[i * 4000:4000 + i * 4000]
        r[:, i, 2] = z[i * 4000:4000 + i * 4000]
        types_n[:, i] = types[i * 4000:4000 + i * 4000]
        atom_idx_n[:, i] = atom_idx[i * 4000:4000 + i * 4000]

    sort_idx = np.argsort(atom_idx_n, axis=0)
    for i in range(r.shape[-1]):
        r[:, :, i] = np.take_along_axis(r[:, :, i], sort_idx, axis=0)

    types_n = np.take_along_axis(types_n, sort_idx, axis=0)

    return r, types_n


def get_dump(filename):
    data = pd.read_table(filename, on_bad_lines="skip", comment="I", skiprows=9, delimiter=" ",
                         names=["atom_nr", "type", "x", "y", "z", "fx", "fy", "fz", "pot_e"]).dropna()
    L = 31.82
    # retrieve unscaled coordinates
    x = np.asarray(data["x"])
    y = np.asarray(data["y"])
    z = np.asarray(data["z"])

    # wrap coordinates
    x = x - np.round(x / L) * L
    y = y - np.round(y / L) * L
    z = z - np.round(z / L) * L
    print(x.shape)

    types = np.asarray(data["type"], dtype=int)
    atom_idx = np.asarray(data["atom_nr"], dtype=int)

    n_frames = int(len(x) // 4000)
    r = np.zeros((4000, n_frames, 3))
    types_n = np.zeros((4000, n_frames), dtype=int)
    atom_idx_n = np.zeros((4000, n_frames), dtype=int)

    for i in range(n_frames):
        r[:, i, 0] = x[i * 4000:4000 + i * 4000]
        r[:, i, 1] = y[i * 4000:4000 + i * 4000]
        r[:, i, 2] = z[i * 4000:4000 + i * 4000]
        types_n[:, i] = types[i * 4000:4000 + i * 4000]
        atom_idx_n[:, i] = atom_idx[i * 4000:4000 + i * 4000]

    sort_idx = np.argsort(atom_idx_n, axis=0)
    for i in range(r.shape[-1]):
        r[:, :, i] = np.take_along_axis(r[:, :, i], sort_idx, axis=0)

    types_n = np.take_along_axis(types_n, sort_idx, axis=0)

    return r, types_n


def cm_movement(x, y, z, types):
    m1 = 26.982
    m2 = 15.999

    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    r[r > 1] -= 1

    masses = np.full_like(types, fill_value=m2)
    masses[types == 1] = m1

    cm_pos = 1 / (np.sum(masses, axis=1)) * np.sum(masses * r, axis=1)
    return cm_pos


def compute_msd(r):
    """
    Computes msd
    :param r: scaled r coordinates, 3d array: (n_atoms, n_frames, n_coordinates = 3).
    :return: msd values for all time intervals possible with the given time frames
            1d array.
    """
    res = np.zeros((r.shape[1] - 1, 2))
    L = 1
    frames = r

    for dt in range(1, frames.shape[1]):  # go through all time intervalls
        for t0 in range(frames.shape[1] - dt):
            dr = frames[:, t0 + dt] - frames[:, t0]
            dr = dr - np.round(dr / L) * L  # pbc
            dr2 = np.sum(dr ** 2, axis=0)
            res[dt - 1, 1] += np.mean(dr2)

        res[dt - 1, 0] = dt
        res[dt - 1, 1] = res[dt - 1, 1] / (frames.shape[1] - dt)
    return res


# r, types = get_data("C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_output\\a_al2o3_4000_1_dump.lammpstrj")
#r, types = get_data("C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_output\\a_al2o3_1.lammpstrj")


def test_msd(r):
    dt = 100  # fs
    # frames_5000T = np.arange(0, 100, 1, dtype=int)
    # frames_3000T = np.arange(300, 400, 1, dtype=int)
    # frames_300T = np.arange(4000, 4100, 1, dtype=int)
    frames_300T = np.arange(250, 260, 1, dtype=int)
    frames_melt = np.arange(0, 50, 1, dtype=int)

    # msd5000T = compute_msd(r[:, frames_5000T])
    # msd3000T = compute_msd(r[:, frames_3000T])
    msd300T = compute_msd(r[:, frames_300T])
    msdmelt = compute_msd(r[:, frames_melt])

    plt.figure()
    plt.title("msd, 300K")
    plt.ylabel("msd")
    plt.xlabel("Time interval [ps]")
    plt.plot(msd300T[:, 0] * 0.1, msd300T[:, 1], "-o")
    plt.show()

    plt.figure()
    plt.plot(msdmelt[:, 0] * 0.1, msdmelt[:, 1], "-o", label="2345K")
    plt.plot(msd300T[:, 0] * 0.1, msd300T[:, 1], "-d", label="300K")
    # plt.plot(msd5000T[:, 0]*0.1, msd5000T[:, 1], label="5000K")
    # plt.plot(msd3000T[:, 0]*0.1, msd3000T[:,1], "--", label="3000K")
    # plt.plot(msd300T[:, 0]*0.1, msd300T[:,1], ":", label="300K")
    plt.title("msd at different temperatures")
    plt.ylabel("msd")
    plt.xlabel("Time interval [ps]")
    plt.legend()
    plt.show()


# test_msd(r)

def make_bonds(r, name_array, box_size, cutoffs):
    """
    Finds atoms connected by chemical bonds and stores the information in matrixes.
    If atom 1 and atom 2 are bonded, then the value at index 1, 2 in the matrix
    will be non-zero.
    Similar bond matrixes will be made to contain bondlengths of the chemical bonds.

    Parameters
    ----------
    x_array : numpy array
        x-coordinates of atoms.
    y_array : numpy array
        y-coordinates of atoms.
    z_array : numpy array
        z-coordinates of atoms.
    name_array : numpy array
        names of atoms. Here: Si or O
    boxx : float
        box size along x-axis.
    boxy :  float
        box size along y-axis.
    boxz :  float
        box size along z-axis.
    cutoffs : list
        Contains the cut-offs which determines if the interatomic distance is
        short enough to be a chemical bonds. [Si-Si, Si-O, O-O]

    Returns
    -------
    index_2 : 2d numpy array
        Number of rows and columns equal to number of atoms.
        If atom 1 and atom 2 are bonded, then the value at index 1, 2 in the matrix
        will be non-zero.
    2d numpy array
        Contains euclidian distances between bonded atoms.
    dist : 2d numpy array
        At index 1, 2 is the vectorized distance between atom 1 and atom 2.

    """
    cut_alal = cutoffs[0]
    cut_alo = cutoffs[1]
    cut_oo = cutoffs[2]

    INDEX = np.arange(0, len(name_array), 1, dtype=int)
    dist = np.zeros((3, len(name_array), len(name_array)))
    index_2 = np.zeros((len(name_array), len(name_array)), dtype=int)

    OO = np.full(len(name_array), 1)
    ALO = np.full(len(name_array), 2)
    ALAL = np.full(len(name_array), 3)

    # bondlengths = []
    for i in range(len(name_array)):
        r1 = r[i]

        r2 = np.copy(r)

        # find distance between the atoms (with periodic boundary):
        dr = r1 - r2 - box_size * np.rint((r1 - r2) / box_size)
        dx, dy, dz = dr[:, 0], dr[:, 1], dr[:, 2]
        dr = np.sqrt(np.sum(dr ** 2, axis=-1))

        dr[i] = np.inf
        values_eq = False
        values_neq = False

        if name_array[i] == 1:
            values_eq = np.logical_and(name_array == name_array[i], dr < cut_alal)
            values_neq = np.logical_and(name_array != name_array[i], dr < cut_alo)

            if np.count_nonzero(values_neq) != 4:  # DEAL WITH THIS!!!!
                # dr = (r[index] - r - box * np.rint((r[index] - r) / box))
                # dr = np.sqrt(np.sum(dr ** 2, axis=2))
                bonds = np.sort(np.copy(dr))[1:5]
                n_idx = np.where(dr == bonds[-1])
                values_neq[n_idx] = True

            btype = ALAL * values_eq + ALO * values_neq
            btype = btype[btype != 0]

        elif name_array[i] == 2:
            values_eq = np.logical_and(name_array == name_array[i], dr < cut_oo)
            values_neq = np.full(values_eq.shape, False)

            btype = OO[values_eq]
        else:
            print('Somthing is wrong with data file')

        val = np.logical_or(values_eq, values_neq)
        index_2[i, INDEX[val]] = btype
        # bondlengths.append(dr[val])
        dist[:, i] = np.array([dx, dy, dz])

    return index_2, dist  # , np.array(bondlengths, dtype=object), dist


def test_make_bonds(r, name_array, box_size, cutoffs):
    frame = -1
    atom_nr1 = 473
    atom_nr2 = 1340
    print(r.shape)
    print("scaled coordinates of atom", atom_nr1)
    print(r[atom_nr1, frame, :])
    print("\nunscaled coordinates of atom", atom_nr1)
    print(r[atom_nr1, frame, :] * box_size)

    print("\n\nscaled coordinates of atom", atom_nr2)
    print(r[atom_nr2, frame, :])
    print("\nunscaled coordinates of atom", atom_nr2)
    print(r[atom_nr2, frame, :] * box_size)

    print("\nunscaled distance between atom", atom_nr1, "and ", atom_nr2)
    print(np.sqrt(np.sum((r[atom_nr2, frame] - r[atom_nr1, frame]) ** 2)) * box_size[0])

    idx, dist = make_bonds(r[:, frame] * box_size[:, None].T, name_array[:, frame], box_size, cutoffs)

    print("dist:")
    dr1 = dist[:, atom_nr1, atom_nr2]
    print(dr1)
    print("r: ", np.sqrt(np.sum(dr1 ** 2, axis=0)))

    bonds = np.sqrt(np.sum((dist ** 2), axis=0)) - np.sum(
        box_size[:, None, None] * np.rint((np.sqrt(dist ** 2)) / box_size[:, None, None]), axis=0)
    print(bonds[atom_nr1, atom_nr2])

    binwidth = 0.1
    plt.figure()
    plt.title("Bondlength distributions")
    plt.hist(bonds[np.triu(idx) == 1], alpha=0.5,
             bins=np.arange(min(bonds[np.triu(idx) == 1]), max(bonds[np.triu(idx) == 1]) + binwidth, binwidth),
             color="darkblue", edgecolor="k", label="Al-Al")
    plt.hist(bonds[np.triu(idx) == 3], alpha=0.5,
             bins=np.arange(min(bonds[np.triu(idx) == 3]), max(bonds[np.triu(idx) == 3]) + binwidth, binwidth),
             color="red", edgecolor="k", label="O-O")
    plt.hist(bonds[np.triu(idx) == 2], alpha=0.5,
             bins=np.arange(min(bonds[np.triu(idx) == 2]), max(bonds[np.triu(idx) == 2]) + binwidth, binwidth),
             color="purple", edgecolor="k", label="Al-O")
    plt.xlabel("Bondlength [Å]")
    plt.legend()
    plt.show()
    print("Al-Al: ", np.count_nonzero(np.triu(idx) == 1))
    print("Al-O: ", np.count_nonzero(np.triu(idx) == 2))
    print("O-O: ", np.count_nonzero(np.triu(idx) == 3))

    half_idx = int(idx.shape[0] // 2)
    bond_count = np.count_nonzero(idx[:half_idx] == 2, axis=0)

    plt.figure()
    plt.title("Number of bonds per atom type")
    plt.hist(bond_count, alpha=0.5, color="darkblue", edgecolor="k", label="Al-O")
    plt.xlabel("Bond count")
    plt.legend()
    plt.show()

    count = []

    for i in range(len(idx)):
        if name_array[i, frame] == 1:
            atomi = np.count_nonzero(idx[i] == 2)
            count.append(atomi)

    print((count[:10]))
    print(np.unique(count, return_counts=True))

    print(np.unique(bond_count, return_counts=True))
    unique, count = np.unique(bond_count, return_counts=True)
    n_Al = 1600

    # Check that numbers are correct
    n_alo_bonds = np.count_nonzero(np.triu(idx) == 2)
    print(n_alo_bonds == np.sum(unique * count))

    for i in range(len(unique)):
        print("Number of bonds: ", unique[i], "\tCount: ", count[i], "\tPercentage: ",
              np.round(count[i] / n_alo_bonds * 100, 2))


L = 31.82


# amorphous:
# test_make_bonds(r, types, box_size=np.array([L, L, L]), cutoffs=[4.5, 2.5, 3.5])
# crystal: did not work
# r_c, types_c = get_data("C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_output\\a_al2o3_4000_1_dump.lammpstrj")
# test_make_bonds(r_c, types_c, box_size=np.array([L, L, L]), cutoffs=[4.5, 2.5, 3.5])

def match_bonds(index):
    """
    Based on the index matrix this function matches bonds of two atoms
    to make chains of three atoms. For calculating angles later
    The different chain configurations have been seperated for
    convienience

    Parameters
    ----------
    index : 2d numpy array
        Number of rows and columns equal to number of atoms.
        If atom 1 and atom 2 are bonded, then the value at index 1, 2 in the matrix
        will be non-zero.

    Returns
    -------
    sisisi_index : 3d numpy array
        If atom 1, 2 and 3, is a Si-Si-Si chain, then at index 1, 2, 3
        there will be a non zero value in the matrix.
    siosi_index : 3d numpy array
        If atom 1, 2 and 3, is a Si-O-Si chain, then at index 1, 2, 3
        there will be a non zero value in the matrix.
    osio_index : 3d numpy array
        If atom 1, 2 and 3, is a O-Si-O chain, then at index 1, 2, 3
        there will be a non zero value in the matrix.

    """
    alalal_index = np.zeros((len(index), len(index), len(index)), dtype=bool)
    aloal_index = np.zeros((len(index), len(index), len(index)), dtype=bool)
    oalo_index = np.zeros((len(index), len(index), len(index)), dtype=bool)

    ALALAL = 7  # Binary 111
    ALOAL = 5  # binary 101
    OALO = 2  # binary (0)10

    ALO = 2
    ALAL = 3

    for i in range(len(index)):
        # types = names[i] #all connected atoms bondtype
        indexes = index[i]  # all connected atoms index

        alal_idx = np.where(indexes == ALAL)[0]  # Finds all Si atoms connected to first Si atom in chain
        alo_idx = np.where(indexes == ALO)[0]  # Finds all O atoms connected to first Si atom in chain
        middle_al = np.where(index[:, i] == ALO)[0]  # Finds all Si atoms connected to first O atom

        for j in range(len(alal_idx)):
            last_al = np.where(index[alal_idx[j]] == ALAL)[
                0]  # finds all Si atoms connected to the middle atom in chain
            last_al = last_al[last_al != i]

            alalal_index[i][alal_idx[j]][last_al] = True

        for j in range(len(alo_idx)):
            last_al = np.where(index.T[alo_idx[j]] == ALO)[
                0]  # Finds all Si atoms connected to the middle atom (O-Si bonds are Si-O bonds reversed -> .T)
            last_al = last_al[last_al != i]

            aloal_index[i][alo_idx[j]][last_al] = True

        for j in range(len(middle_al)):
            last_o = np.where(index[middle_al[j]] == ALO)[0]
            last_o = last_o[last_o != i]

            oalo_index[i][middle_al[j]][last_o] = True

    return alalal_index, aloal_index, oalo_index


def calc_angle(index_matrix, dr):
    """
    Uses the cosine law to calculate the angle between 3 atoms using
    bond lengths.

    Parameters
    ----------
    index_matrix : 3d numpy array of ints and 0s.
        If index_matrix[i, j, k] != 0 there is a bond between atom i, j and j, k.
        The value at this place is different for different bonds.
    dr : 3d numpy array of floats.
        At dr[i, j] we find the vector between atom i and j.

    Returns
    -------
    angle : 1d numpy array, floats
        Array of angles.

    """
    indexes = np.array(np.where(index_matrix != 0)).T

    index1 = indexes[:, :2]
    index2 = indexes[:, 1:]
    a = np.sqrt(np.sum(dr[:, index1[:, 0], index1[:, 1]] ** 2, axis=0))
    b = np.sqrt(np.sum(dr[:, index2[:, 0], index2[:, 1]] ** 2, axis=0))
    c = np.sqrt(np.sum(dr[:, index1[:, 0], index2[:, 1]] ** 2, axis=0))

    angle = np.rad2deg(np.arccos((a ** 2 + b ** 2 - c ** 2) / (2 * a * b)))  # Cosine law ##
    return angle


def get_angles(dr, name_array, cutoffs, bz=1):
    n = int(180//bz)+1
    alalal = np.zeros(n)
    aloal = np.zeros(n)
    oalo = np.zeros(n)
    ooo = np.zeros(n)
    angle_ticks = np.linspace(0, 180, n)

    for i in range(dr.shape[0]):
        typei = name_array[i]

        for j in range(dr.shape[0]):
            rij = np.sqrt(np.sum((dr[i] - dr[j])**2))  # distance between atom i and j
            if i == j:
                continue  # skip this iteration
            typej = name_array[j]

            if rij <= max(cutoffs):  # distance must be within cut off radius
                for k in range(dr.shape[0]):
                    if i == k or j == k:
                        continue  # skip this iteration

                    rkj = np.sqrt(np.sum((dr[j] - dr[k])**2))  # distance between atom j and k
                    typek = name_array[k]

                    if rkj <= max(cutoffs):  # distance must be within cut off radius
                        rik = np.sqrt(np.sum((dr[i] - dr[k])**2))

                        angle = np.rad2deg(np.arccos((rij ** 2 + rkj ** 2 - rik ** 2) / (2 * rij * rkj)))  # Cosine law

                        if np.isnan(angle):  # Skip if there is no angle
                            continue

                        idx = int(np.rint(angle/bz))

                        if rij <= cutoffs[1] and rkj <= cutoffs[1]:
                            if typei == 1 and typej == 2 and typek == 1:  # AlOAl
                                aloal[idx] += 1
                            elif typei == 2 and typej == 1 and typek == 2:  # OAlO
                                oalo[idx] += 1
                        elif (typei == 1 and typej == 1 and typek == 1) and (
                                rij <= cutoffs[0] and rkj <= cutoffs[0]):  # AlAlAl
                            alalal[idx] += 1
                        elif (typei == 2 and typej == 2 and typek == 2) and (
                                rij <= cutoffs[2] and rkj <= cutoffs[2]):  # OOO
                            ooo[idx] += 1

    return alalal, aloal, oalo, ooo, angle_ticks

def test_angles2():
    #f100 = "C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_dump\\a_al2o3_1.dump"
    f100 = "../lammps/output/a_al2o3_1.dump"
    f10 = "../lammps/output/a_al2o3_23.dump"
    f200 = "../lammps/output/a_al2o3_55.dump"
    dr100, name_array100 = get_dump(f100)
    dr10, name_array10 = get_dump(f10)
    dr200, name_array200 = get_dump(f200)

    cutoffs = [4, 2.5, 3.4]
    frame = -1

    alalal100, aloal100, oalo100, ooo100, angle_ticks = get_angles(dr100[:, frame], name_array100[:, frame], cutoffs, bz=5)
    alalal10, aloal10, oalo10, ooo10, angle_ticks = get_angles(dr10[:, frame], name_array10[:, frame], cutoffs, bz=5)
    alalal200, aloal200, oalo200, ooo200, angle_ticks = get_angles(dr200[:, frame], name_array200[:, frame], cutoffs, bz=5)



    plt.figure()
    plt.title("Al-Al-Al adf")
    plt.plot(angle_ticks, alalal100, "-o", label="100K/ps")
    plt.plot(angle_ticks, alalal10, "-d", label="10K/ps")
    plt.plot(angle_ticks, alalal200, "-s", label="200K/ps")
    plt.legend()
    plt.savefig("alalal_adf.png")
    #plt.show()

    plt.figure()
    plt.title("Al-O-Al adf")
    plt.plot(angle_ticks, aloal100, "-o", label="100K/ps")
    plt.plot(angle_ticks, aloal10, "-d", label="10K/ps")
    plt.plot(angle_ticks, aloal200, "-s", label="200K/ps")
    plt.legend()
    plt.savefig("aloal_adf.png")
    #plt.show()

    plt.figure()
    plt.title("O-Al-O adf")
    plt.plot(angle_ticks, oalo100, "-o", label="100K/ps")
    plt.plot(angle_ticks, oalo10, "-d", label="10K/ps")
    plt.plot(angle_ticks, oalo200, "-s", label="200K/ps")
    plt.legend()
    plt.savefig("oalo_adf.png")
    #plt.show()

    plt.figure()
    plt.title("O-O-O adf")
    plt.plot(angle_ticks, ooo100, "-o", label="100K/ps")
    plt.plot(angle_ticks, ooo10, "-d", label="10K/ps")
    plt.plot(angle_ticks, ooo200, "-s", label="200K/ps")
    plt.legend()
    plt.savefig("ooo_adf.png")
    #plt.show()

#test_angles2()

def test_angles():
    r, name_array = get_dump("C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_dump\\a_al2o3_1.dump")

    L = 31.82
    box_size = np.array([L, L, L])

    cutoffs = [4.5, 2.5, 3.5]
    frame = -1
    idx, dist = make_bonds(r[:, frame] * box_size[:, None].T, name_array[:, frame], box_size, cutoffs)
    alalal_index, aloal_index, oalo_index = match_bonds(idx)
    aloal_angle = calc_angle(aloal_index, dist)
    oalo_angle = calc_angle(oalo_index, dist)

    binwidth = 0.1
    plt.figure()
    plt.title("Bondangle distributions")
    plt.hist(aloal_angle, alpha=0.5, bins=np.arange(min(aloal_angle), max(aloal_angle) + binwidth, binwidth),
             color="darkblue", edgecolor="k", label="Al-Al")
    plt.hist(oalo_angle, alpha=0.5, bins=np.arange(min(oalo_angle), max(oalo_angle) + binwidth, binwidth), color="red",
             edgecolor="k", label="O-O")
    # plt.hist(bonds[np.triu(idx) == 2], alpha=0.5, bins=np.arange(min(bonds[np.triu(idx) == 2]), max(bonds[np.triu(idx) == 2]) + binwidth, binwidth), color="purple", edgecolor="k", label="Al-O")
    plt.xlabel("Bondangle [degree]")
    plt.legend()
    plt.show()


#test_angles()


def plot_cm():
    cm_pos = cm_movement(x, y, z, types)
    cm_time = np.linspace(0, 410, len(cm_pos))

    plt.figure()
    plt.title("Center of mass over time")
    plt.plot(cm_time, cm_pos, label="CM position")
    plt.plot(cm_time, np.full_like(cm_time, np.mean(cm_pos)), label="Average CM position")
    plt.xlabel("Time [ps]")
    plt.ylabel("Position/box-dimension")
    plt.legend()
    plt.grid()
    plt.show()


def rdf(filename):
    # filename = "C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_output\\alo4000_2_melt.rdf"
    rdf = pd.read_table(filename, on_bad_lines="skip", comment="#", skiprows=4, delimiter=" ",
                        names=[1, 2, 3, 4, 5, 6, 7, 8]).dropna()  # .reset_index(drop=True)

    dist = np.asarray(rdf[2])
    rdf1 = np.asarray(rdf[3])  # Al Al
    rdf2 = np.asarray(rdf[4])  # other thing
    rdf3 = np.asarray(rdf[5])  # Al O
    rdf4 = np.asarray(rdf[6])  # other thing
    rdf5 = np.asarray(rdf[7])  # O O
    rdf6 = np.asarray(rdf[8])  # other thing

    from matplotlib.animation import FuncAnimation
    from matplotlib import animation

    def animate_rdf():
        # data which the line will
        # contain (x, y)
        def init():
            line.set_data([], [])
            return line,

        def animate(i):
            x = np.array(dist[:1000])
            y = np.array(rdf5[i * 1000:1000 + i * 1000])
            if len(y) == 0:
                print(i)
            line.set_data(x, y)
            title.set_text("T = " + str(3000 - 10 * i))
            return line, title,

        n_frames = 270

        # initializing a figure in
        # which the graph will be plotted
        fig = plt.figure()

        # marking the x-axis and y-axis
        axis = plt.axes(xlim=(0, 10), ylim=(0, 6))
        plt.title("O-O rdf")
        plt.xlabel("Distance [Å]")
        plt.ylabel("g(r)")
        ax = plt.gca()

        # initializing a line variable
        line, = axis.plot([], [], lw=2)
        title = ax.text(0.9, 0.9, "T=3000K", bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5},
                        transform=ax.transAxes, ha="center")

        # anim = FuncAnimation(fig, animate, init_func=init,
        #                     frames=n_frames, interval=100, blit=True)

        # writervid = animation.PillowWriter(fps=10)
        # anim.save("OO.gif", writer=writervid)

    plt.figure()
    plt.title("RDF")
    plt.plot(dist[:1000], rdf1[-1000:], label="Al-Al", color="darkblue")
    plt.plot(dist[:1000], rdf3[-1000:], label="Al-O", color="purple")
    plt.plot(dist[:1000], rdf5[-1000:], label="O-O", color="red")
    plt.legend()
    plt.xlabel("Distance [Å]")
    plt.ylabel("g(r)")
    plt.show()

    plt.figure()
    plt.title("RDF Al-O, T = 300K")
    plt.plot(dist[:1000], rdf3[-1000:])
    plt.xlabel("Distance [Å]")
    plt.ylabel("g(r)")
    plt.show()

    plt.figure()
    plt.title("RDF O-O, T = 300K")
    plt.plot(dist[:1000], rdf5[-1000:])
    plt.xlabel("Distance [Å]")
    plt.ylabel("g(r)")
    plt.show()


def rdf_cooling(filename, c_rate="100K/ps"):
    rdf = pd.read_table(filename, on_bad_lines="skip", comment="#", skiprows=4, delimiter=" ",
                        names=[1, 2, 3, 4, 5, 6, 7, 8]).dropna()  # .reset_index(drop=True)
    dist = np.asarray(rdf[2]).reshape(-1, 1000)
    rdf1 = np.asarray(rdf[3]).reshape(-1, 1000)  # Al Al
    rdf3 = np.asarray(rdf[5]).reshape(-1, 1000)  # Al O
    rdf5 = np.asarray(rdf[7]).reshape(-1, 1000)  # O O

    # idx = dist.argsort()
    plt.figure()
    plt.title("RDF, cooling rate = " + c_rate)
    plt.plot(np.mean(dist, axis=0), np.mean(rdf1, axis=0), label="Al-Al", color="darkblue")
    plt.plot(np.mean(dist, axis=0), np.mean(rdf3, axis=0), label="Al-O", color="purple")
    plt.plot(np.mean(dist, axis=0), np.mean(rdf5, axis=0), label="O-O", color="red")
    plt.xlabel("Distance [Å]")
    plt.legend()
    plt.ylabel("g(r)")
    plt.show()

    return np.mean(dist, axis=0), np.mean(rdf1, axis=0)


# rdf_cooling("C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_output\\a_al2o3_100_Kps.rdf")
# dist, rdf200 = rdf_cooling("C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_output\\a_al2o3_200Kps.rdf", c_rate="200K/ps")
# dist, rdf100 = rdf_cooling("C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_output\\a_al2o3_100_Kps.rdf", c_rate="100K/ps")
# dist, rdf10 = rdf_cooling("C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_output\\a_al2o3_10Kps.rdf", c_rate="10K/ps")


# plt.figure()
# plt.title("Al-Al rdf for different cooling rates")
# plt.plot(dist, rdf200, label="200K/ps")
# plt.plot(dist, rdf10, label="10K/ps")
# plt.xlabel("Distance [Å]")
# plt.ylabel("g(r)")
# plt.xlim(0, 2)
# plt.ylim(-0.01, 0.01)
# plt.grid()
# plt.legend()
# plt.show()

"""
rdf = get_rdf("C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_output\\alo4000.rdf")
dist = np.asarray(rdf[2])
rdf1 = np.asarray(rdf[3])
rdf2 = np.asarray(rdf[4])

plt.figure()
plt.plot(dist[:1000], rdf1[:1000])
plt.show()

plt.figure()
plt.title("RDF, 4000 atoms")
plt.plot(dist[:1000], np.mean(np.reshape(rdf1, (-1, 1000, 1)), axis=0))
plt.xlabel("Distance [Å]")
plt.ylabel("g(r)")
plt.show()


filename = "C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_output\\alo4000_1.rdf"
rdf = pd.read_table(filename, on_bad_lines="skip", comment="#", skiprows=4, delimiter=" ",
                    names=[1, 2, 3, 4, 5, 6, 7, 8]).dropna()  # .reset_index(drop=True)

dist = np.asarray(rdf[2])
rdf1 = np.asarray(rdf[3]) #Al Al
rdf2 = np.asarray(rdf[4]) #other thing
rdf3 = np.asarray(rdf[5]) #Al O
rdf4 = np.asarray(rdf[6]) #other thing
rdf5 = np.asarray(rdf[7]) #O O
rdf6 = np.asarray(rdf[8]) #other thing

plt.figure()
plt.title("RDF Al-Al, 4000 atoms")
plt.plot(dist[:1000], np.mean(np.reshape(rdf1, (-1, 1000, 1)), axis=0))
plt.xlabel("Distance [Å]")
plt.ylabel("g(r)")
plt.show()

plt.figure()
plt.title("RDF Al-O, 4000 atoms")
plt.plot(dist[:1000], np.mean(np.reshape(rdf3, (-1, 1000, 1)), axis=0))
plt.xlabel("Distance [Å]")
plt.ylabel("g(r)")
plt.show()


plt.figure()
plt.title("RDF O-O, 4000 atoms")
plt.plot(dist[:1000], np.mean(np.reshape(rdf5, (-1, 1000, 1)), axis=0))
plt.xlabel("Distance [Å]")
plt.ylabel("g(r)")
plt.show()
"""
