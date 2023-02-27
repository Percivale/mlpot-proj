import numpy as np

def read_dump(path: str):
    """
    :param path: path to dump file
    :param n_atoms: number of atoms in simulation
    :return:
    """
    # Maybe make the below line general, and take usecols as a function argument??
    xyz = np.genfromtxt(path, dtype=float, comments="ITEM", usecols=[2, 3, 4],
                        invalid_raise=False, skip_header=9)
    print(xyz[0])

    idx_type = np.genfromtxt(path, dtype=int, comments="ITEM", usecols=[0, 1],
                        invalid_raise=False, skip_header=9)

    print(xyz.shape, xyz[0, 0, :])
    print(idx_type.shape, idx_type[0, 0, :])

    boxsize = np.genfromtxt(path, dtype=float, comments="ITEM", skip_header=3, max_rows=3)
    n_atoms = np.loadtxt(path, dtype=int, comments="ITEM", skiprows=3, max_rows=1)[0]

    xyz = np.reshape(xyz, (-1, n_atoms, 3))
    idx_type = np.reshape(idx_type, (-1, n_atoms, 2))

    # TO DO:
    # get forces and energies??

    return xyz, idx_type, boxsize, n_atoms


def format_files():
    # loop through files in directory
    # read dump files
    # reformat relevant data
    # save in new format
    return 0


def test():
    path = "C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_output\\a_al2o3_0_dump.lammpstrj"
    xyz, idx_type, boxsize, n_atoms = read_dump(path)

    print("TEST")
    print(xyz.shape, xyz[0, 0, :])
    print(idx_type.shape, idx_type[0, 0, :])
    print(boxsize)
    print(n_atoms)


test()
