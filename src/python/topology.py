import numpy as np
import matplotlib.pyplot as plt

def make_top(n1:int, n2:int, q1:float, q2:float, n_atoms:int, sep:float):
    boxsize = np.round(n_atoms**(1/3)*sep+1)
    l = np.arange(0.5, boxsize, sep)
    xyz = np.vstack(np.meshgrid(l, l, l)).reshape(3,-1).T

    idx = np.arange(1, n_atoms+1, 1, dtype=int)

    type = np.concatenate((np.ones(n1, dtype=int), np.ones(n2, dtype=int)*2))
    q = np.concatenate((np.ones(n1, dtype=float)*q1, np.ones(n2, dtype=float)*q2))

    for i in range(n_atoms):
       print(idx[i], type[i], q[i], xyz[i, 0], xyz[i, 1], xyz[i, 2])


def make_alo(n1:int, n2:int, q1:float, q2:float, a:float, c:float, angle=60):
    #make unit cell:
    xycell = np.array([[0, 0, 0],
                     [0, 1, 0],
                     [np.sqrt(3)/2, -1/2, 0],
                     [np.sqrt(3)/2, 1/2, 0]])
    xycell[:, :-1] *= a
    xycellc = np.copy(xycell)
    xycellc[:, -1] += c
    xyzcell = np.concatenate((xycell, xycellc), axis=0)

    #make one hexagon
    xyz_hexa = hexa(xyzcell)
    print(xyz_hexa.shape)

    #make many hexagons
    n_cells = int(np.rint((n1+n2)/14))
    xyz = np.copy(xyz_hexa)
    print(n_cells)
    print(xyz_hexa.shape[0]*n_cells)

    types = np.zeros(xyz_hexa.shape[0]*n_cells+1, dtype=str)
    types[:] = "red"
    type_idx = np.array([1, 3, 5, 7, 9, 11, 13, 15, 18, 19, 22, 23])
    types[type_idx] = "blue"
    print(types)

    plt.figure()
    for i in range(xyz_hexa.shape[0]):
        plt.plot(xyz_hexa[i, 0], xyz_hexa[i, 1], "o", color=types[i])
    plt.show()

    for i in range(n_cells):
        dir = i%3
        temp = np.copy(xyz_hexa)
        temp[:, dir] += (np.max(temp[:, dir]) - np.min(temp[:, dir]))
        xyz = np.concatenate((xyz, temp), axis=0)

        types[type_idx+n_cells*i] = "blue"

    xyz = np.unique(xyz, axis=0)
    plt.figure()
    for i in range(xyz.shape[0]):
        plt.plot(xyz[i, 0], xyz[i, 1], "o", color=types[i])
    plt.show()


    return 0

def hexa(xyz):
    temp = np.copy(xyz)
    left = np.zeros_like(xyz)
    left[:, 0] -= temp[:, 0]
    left[:, 1:] += temp[:, 1:]
    down = np.zeros_like(xyz)
    down[:, 1] -= 2*temp[:, 1]
    down += temp
    xyz = np.concatenate((left, down, xyz), axis=0)

    return np.unique(xyz, axis=0)


def testalo():
    n1 = 20
    n2 = 30
    a = 4.785
    c = 12.991
    q1 = 1.4175
    q2 = -0.9450
    make_alo(n1, n2, q1, q2, a, c)

testalo()


def test():
    n1 = 200
    n2 = 300
    n_atoms = 500
    sep = 1.5
    q1 = 1.4175
    q2 = -0.9450
    make_top(n1, n2, q1, q2, n_atoms, sep)

#test()
def rep_pot():
    r = np.linspace(0.01, 1.0, 100)
    pot = 15/r**24
    f = 15*24/r**25

    for i in range(len(r)):
        print(i+1, r[i], pot[i], f[i])


