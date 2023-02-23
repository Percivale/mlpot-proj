import numpy as np

def read_dump(path: str, n_atoms: int):
    """
    :param path: path to dump file
    :param n_atoms: number of atoms in simulation
    :return:
    """
    # Maybe make the below line general, and take usecols as a function argument??
    xyz = np.genfromtxt(path, dtype=float, comments="ITEM", usecols=[2, 3, 4],
                        invalid_raise=False, skip_header=9)
    print("TEST")
    print(xyz.shape)
    print(xyz[0:3])
    xyz = np.reshape(xyz, (-1, 5000, 3))
    print(xyz[0, 0:3, :])
    print(xyz.shape)

    # TO DO:
    # get box size
    # get n_atoms
    # rescale xyz ?
    # get index and names data
    # get forces and velocities?

    return xyz


def format_files():
    # loop through files in directory
    # read dump files
    # reformat relevant data
    # save in new format
    return 0


def test():
    path = "C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_output\\a_al2o3_0_dump.lammpstrj"
    read_dump(path, 5000)


test()
