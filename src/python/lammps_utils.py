import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

class ConvertUnits:
    def __init__(self):
        self.kcal_eV = 3.82929406e-23
        self.mol = 6.02214076e23
        self.kcal_kJ = 4.184
        self.u = 1.66e-27
        self.Å = 1e-10

    def eV_to_kcal_mol(self, num):
        num = num * self.kcal_eV * self.mol
        return num

    def kcal_mol_to_eV(self, num):
        num = num / self.kcal_eV / self.mol
        return num

    def kcal_to_kJ(self, num):
        num = num * self.kcal_kJ
        return num

    def kJ_to_kcal(self, num):
        num = num / self.kcal_kJ
        return num

    def gram_cm3_to_u_Å3(self, num):
        num = num * 1000 * self.mol / (self.Å ^ 3)
        return num

CU = ConvertUnits()
print(CU.kcal_mol_to_eV(-741073.45))
print(CU.kcal_mol_to_eV(23.069))
print(CU.kcal_mol_to_eV(-762301.51))

def test_CU():
    CU = ConvertUnits()
    print(CU.kcal_mol_to_eV(1))
    print(CU.eV_to_kcal_mol(1 / 55.26349406) * 10000)
    print(CU.eV_to_kcal_mol(0.00294848))
    print(CU.eV_to_kcal_mol(0.00745792))
    print(CU.eV_to_kcal_mol(0.01196736))
    print(CU.eV_to_kcal_mol(14.561856))
    print(CU.eV_to_kcal_mol(35.836056))
    print(CU.eV_to_kcal_mol(88.190881))
    print(CU.kcal_mol_to_eV(324.0152636))
    print(CU.kcal_mol_to_eV(797.3661878))
    print(CU.kcal_mol_to_eV(1962.2311319))

#test_CU()


def read_dump(filename: str, n_atoms=4000):
    data = pd.read_table(filename, on_bad_lines="skip", comment="I", skiprows=9, delimiter=" ",
                         names=["atom_nr", "type", "x", "y", "z", "fx", "fy", "fz", "pe"]).dropna()  # .reset_index(drop=True)
    x = np.asarray(data["x"])
    y = np.asarray(data["y"])
    z = np.asarray(data["z"])
    #vx = np.asarray(data["vx"])
    #vy = np.asarray(data["vy"])
    #vz = np.asarray(data["vz"])
    fx = np.asarray(data["fx"])
    fy = np.asarray(data["fy"])
    fz = np.asarray(data["fz"])
    pe = np.asarray(data["pe"])
    types = np.asarray(data["type"], dtype=int)
    atom_idx = np.asarray(data["atom_nr"], dtype=int)

    n_frames = int(len(x) // n_atoms)
    r = np.zeros((n_atoms, n_frames, 3))
    f = np.zeros((n_atoms, n_frames, 3))
    #v = np.zeros((n_atoms, n_frames, 3))
    types_n = np.zeros((n_atoms, n_frames), dtype=int)
    atom_idx_n = np.zeros((n_atoms, n_frames), dtype=int)

    for i in range(n_frames):
        r[:, i, 0] = x[i * n_atoms:n_atoms + i * n_atoms]
        r[:, i, 1] = y[i * n_atoms:n_atoms + i * n_atoms]
        r[:, i, 2] = z[i * n_atoms:n_atoms + i * n_atoms]

        #v[:, i, 0] = vx[i * n_atoms:n_atoms + i * n_atoms]
        #v[:, i, 1] = vy[i * n_atoms:n_atoms + i * n_atoms]
        #v[:, i, 2] = vz[i * n_atoms:n_atoms + i * n_atoms]

        f[:, i, 0] = fx[i * n_atoms:n_atoms + i * n_atoms]
        f[:, i, 1] = fy[i * n_atoms:n_atoms + i * n_atoms]
        f[:, i, 2] = fz[i * n_atoms:n_atoms + i * n_atoms]

        types_n[:, i] = types[i * n_atoms:n_atoms + i * n_atoms]
        atom_idx_n[:, i] = atom_idx[i * n_atoms:n_atoms + i * n_atoms]

    sort_idx = np.argsort(atom_idx_n, axis=0)
    for i in range(r.shape[-1]):
        r[:, :, i] = np.take_along_axis(r[:, :, i], sort_idx, axis=0)
        f[:, :, i] = np.take_along_axis(f[:, :, i], sort_idx, axis=0)
        #v[:, :, i] = np.take_along_axis(v[:, :, i], sort_idx, axis=0)
    types_n = np.take_along_axis(types_n, sort_idx, axis=0)

    return r, f, types_n[:, 0], pe


def get_num(lis):
    for num in lis:
        try:
            yield float(num)
        except ValueError:
            pass


def read_log(filename: str):
    file = open(filename, "r")
    cols = ["Step", "Temp", "E_pair", "E_mol", "TotEng", "Press"]
    file_list = file.readlines()
    data = []

    for line in file_list:
        line_list = line.split()
        if len(line_list) == len(cols):
            data_line = list(get_num(line_list))
            if len(data_line) == len(cols):
                data.append(data_line)
    file.close()
    return np.unique(np.asarray(data), axis=0)


def format_logdir(data_dir, dest_dir):
    liquid_frames = [7]
    solid_frames = [-2]
    melt_100Kps = list(np.arange(9, 29, 1, dtype=int))
    melt_10Kps = list(np.arange(9, 213, 1, dtype=int))
    melt_200Kps = list(np.arange(9, 20, 1, dtype=int))

    frames100 = liquid_frames + melt_100Kps + solid_frames
    frames200 = liquid_frames + melt_200Kps + solid_frames
    frames10 = liquid_frames + melt_10Kps + solid_frames
    frames100h = liquid_frames + list(np.arange(9, 32, 1, dtype=int)) + solid_frames
    #frames200h = liquid_frames + list(np.arange(9, 20, 1, dtype=int)) + solid_frames
    frames10h = liquid_frames + list(np.arange(9, 243, 1, dtype=int)) + solid_frames
    CU = ConvertUnits()

    for filename in os.listdir(data_dir):
        f = os.path.join(data_dir, filename)
        if filename[-3:] == "log":
            data = read_log(f)
            if len(data) - 10 == len(frames10):
                energies = data[frames10, 4]
            elif len(data) - 10 == len(frames100):
                energies = data[frames100, 4]
            elif len(data) - 10 == len(frames100h):
                energies = data[frames100h, 4]
            elif len(data) - 10 == len(frames10h):
                energies = data[frames10h, 4]
            #elif len(data) - 10 == len(frames200h):
            #    energies = data[frames200h, 4]
            else:
                energies = data[frames200, 4]
            filename_dir = "/" + filename[:-4] + "/set.000"
            os.makedirs(dest_dir + filename_dir, exist_ok=True)
            dest_path = os.path.join(dest_dir + filename_dir, "energy.npy")
            np.save(file=dest_path, arr=CU.kcal_mol_to_eV(energies))



def format_system(fname):
    xyz, fxyz, types_n, pe = read_dump(fname)
    n_frames = xyz.shape[1]
    n_atoms = xyz.shape[0]

    xyz_raw = np.zeros(shape=(n_frames, n_atoms*3), dtype=float)
    forces_raw = np.zeros(shape=(n_frames, n_atoms*3), dtype=float)
    types_raw = types_n - 1

    for frame_idx in range(n_frames):
        for atom_i in range(n_atoms):
            xyz_raw[frame_idx, atom_i*3:atom_i*3+3] = xyz[atom_i, frame_idx]
            forces_raw[frame_idx, atom_i*3:atom_i*3+3] = fxyz[atom_i, frame_idx]

    return xyz_raw, forces_raw, types_raw


def format_files(data_dir, dest_dir):
    CU = ConvertUnits()
    box_start = 0
    box_end = 31.82

    box_line = np.array([box_end, box_start, box_start, box_start, box_end, box_start, box_start, box_start, box_end])
    for filename in os.listdir(data_dir):
        f = os.path.join(data_dir, filename)
        print("accessing file, ", filename)

        coords, forces, types = format_system(f)
        box = np.stack([box_line]*coords.shape[0], axis=0)


        filename_dir = "/" + filename[2:-5] + "/set.000"

        path_types = os.path.join(dest_dir + "/" + filename[2:-5], "type.raw")
        path_coords = os.path.join(dest_dir + filename_dir, "coord.npy")
        path_forces = os.path.join(dest_dir + filename_dir, "force.npy")
        path_box = os.path.join(dest_dir + filename_dir, "box.npy")

        print("saving forces to: ", path_forces)
        np.save(file=path_forces, arr=CU.kcal_mol_to_eV(forces[4:]))
        print("saving coordinates to: ", path_coords)
        np.save(file=path_coords, arr=coords[4:])
        print("saving types to: ", path_types)
        np.savetxt(path_types, types, fmt='%s', delimiter="\n")
        print("saving box to: ", path_box)
        np.save(file=path_box, arr=box[4:])

#fname = "..\\lammps_script\\a_al2o3_dump\\a_al2o3_5.dump"

#coords, forces, types = format_system(fname)
#print(coords[4:].shape)
#print(forces[4:].shape)
#print(types.shape)

#data = "../deepmd_data/dumps"
#dest = "../deepmd_data/data_set"
#format_files(data, dest)


#dest = "../deepmd_data/data_set"
#data = "../deepmd_data/logs"
#format_logdir(data_dir=".\\log_test", dest_dir=".\\log_test")
#format_logdir(data_dir=data, dest_dir=dest)


# format_log(data_1, melt_frames=melt_100Kps, file_dest="energies_100Kps.npy")  # this looks good
# format_log(data_2, melt_10Kps, file_dest="energies_10Kps.npy")
# format_log(data_3, melt_200Kps, file_dest="energies_200Kps.npy")

def get_adf(filename: str, cr: float):
    adf = pd.read_table(filename, on_bad_lines="skip", comment="#", skiprows=4, delimiter=" ",
                        names=[1, 2, 3, 4, 5, 6]).dropna()  # .reset_index(drop=True)
    print(adf.head(5))
    angle = np.asarray(adf[2]).reshape((-1, 32))
    adf1 = np.asarray(adf[3]).reshape((-1, 32))  # Al O Al
    adf2 = np.asarray(adf[4]).reshape((-1, 32))  # other thing
    adf3 = np.asarray(adf[5]).reshape((-1, 32))  # Al O
    adf4 = np.asarray(adf[6]).reshape((-1, 32))  # other thing

    print(angle[-1])
    print(adf1[-1])
    plt.figure()
    plt.title("Angular distribution function, cr = " + str(cr) + "K/ps")
    plt.plot(np.mean(angle, axis=0), np.mean(adf1, axis=0), label="Al-O-Al")
    # plt.plot(np.mean(angle, axis=0), np.mean(adf2, axis=0), label="idk2")
    plt.plot(np.mean(angle, axis=0), np.mean(adf3, axis=0), label="O-Al-O")
    # plt.plot(np.mean(angle, axis=0), np.mean(adf4, axis=0), label="idk4")
    plt.legend()
    plt.show()

def plot_model(path):
    data = np.genfromtxt(path, names=True)
    for name in data.dtype.names[1:-1]:
        plt.plot(data['step'], data[name], label=name)
    plt.legend()
    plt.xlabel('Step')
    plt.ylabel('Loss')
    plt.xscale('symlog')
    plt.yscale('log')
    plt.grid()
    plt.show()

path = "deepmd_data/lcurve.out"
#plot_model(path)
# get_adf(filename="C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_output\\a_al2o3_1.adf", cr=100)
# get_adf(filename="C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_output\\a_al2o3_23.adf", cr=10)
# get_adf(filename="C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_output\\a_al2o3_45.adf", cr=200)


