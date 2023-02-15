import numpy as np

def read_lammpstrj(path):
    f = open(path, 'r')

    get_val = False
    info = np.zeros(5)
    idx = 0
    for line in f:
        if get_val:
            try:
                val = int(line)
            except ValueError:
                val = float(line)
        if line[:4] == "ITEM":
            get_val = True


    f.close()