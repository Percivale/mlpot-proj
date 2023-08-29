kcal_eV = 3.82929406e-23
mol = 6.02214076e23
kcal_kJ = 4.184
u = 1.66e-27
Å = 1e-10


def eV_to_kcal_mol(num):
    num = num * kcal_eV * mol
    return num

def kcal_mol_to_eV(num):
    num = num / kcal_eV / mol
    return num

def kcal_to_kJ(num):
    num = num * kcal_kJ
    return num

def kJ_to_kcal(num):
    num = num / kcal_kJ
    return num

def gram_cm3_to_u_Å3(num):
    num = num * 1000 * mol / (Å ^ 3)
    return num

def tests():
    print(kcal_mol_to_eV(1))
    print(eV_to_kcal_mol(1 / 55.26349406) * 10000)
    print(eV_to_kcal_mol(0.00294848))
    print(eV_to_kcal_mol(0.00745792))
    print(eV_to_kcal_mol(0.01196736))
    print(eV_to_kcal_mol(14.561856))
    print(eV_to_kcal_mol(35.836056))
    print(eV_to_kcal_mol(88.190881))
    print(kcal_mol_to_eV(324.0152636))
    print(kcal_mol_to_eV(797.3661878))
    print(kcal_mol_to_eV(1962.2311319))
    print(kcal_mol_to_eV(23.069))

