import dpdata
import numpy as np
import os
filename = "C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_dump\\a_al2o3_3.dump"

def coords(filename):
    frames100 = list(np.arange(4, 26, 1, dtype=int))
    frames200 = list(np.arange(4, 210, 1, dtype=int))
    frames10 = list(np.arange(4, 17, 1, dtype=int))

    data = dpdata.System(filename, unwrap=True).sub_system(frames100)
    #data.to("deepmd/npy", "deepmd_data\\a_al2o3_3", set_size=data.get_nframes())
    data.to("deepmd/raw", "deepmd_data\\a_al2o3_3")

coords("..\\lammps_script\\a_al2o3_dump\\a_al2o3_1.dump")
#data = np.load("deepmd_data\\a_al2o3_1\\set.000\\coord.npy")

def get_training(data_dir:str, dest_dir:str, frames=-1):
    #specify frames here.... dobbel check that it is right!!! do it
    frames_melt = [0, 11]
    frames_cooling = [56, -11] #check that these are right!
    frames_solid = [-11, -1]

    for filename in os.listdir(data_dir):
        f = os.path.join(data_dir, filename)
        if os.path.isfile(f):
            data = dpdata.System(filename, unwrap=True)
            data.to("deepmd/npy", os.path.join(dest_dir, filename[:-5]), set_size=data.get_nframes(), frame_idx=frames)

    return 0

#get_training("C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\lammps_script\\a_al2o3_dump", "deepmd_data")

energies = np.load("C:\\Users\\kajah\\git_repo\\mlpot-proj\\src\\data\\energy.npy")

print(energies.shape)
print(energies[:5])
