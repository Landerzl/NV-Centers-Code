import numpy as np
from pyscf import gto, scf
from pyscf.prop import zfs
from pyscf.tools import cubegen

# 1. SISTEMA (NV- Pasivado)
nv_pasivado = '''
N    0.000000    0.000000    0.000000
C    1.020000    1.020000    1.020000
C   -1.020000   -1.020000    1.020000
C    1.020000   -1.020000   -1.020000
H    1.400000    2.000000    0.800000
H    1.800000    0.600000    1.600000
H   -1.400000   -2.000000    0.800000
H   -1.800000   -0.600000    1.600000
H    1.400000   -2.000000   -0.800000
H    1.800000   -0.600000   -1.600000
'''

mol = gto.M(atom=nv_pasivado, basis='def2-svp', charge=-1, spin=2, verbose=4)

# 2. CAMPO ELÉCTRICO 
E_x, E_y, E_z = 0.002, 0.0, 0.0  

# 3. CÁLCULO SCF
mf = scf.UKS(mol)
mf.xc = 'b3lyp'
mf.level_shift = 0.2

hcore_orig = mf.get_hcore
def get_hcore_con_campo(mol=mol):
    h = hcore_orig(mol)
    mol.set_common_orig([0, 0, 0])
    dipole_ints = mol.intor_symmetric('int1e_r', comp=3)
    return h + (E_x * dipole_ints[0] + E_y * dipole_ints[1] + E_z * dipole_ints[2])

mf.get_hcore = get_hcore_con_campo
mf.kernel()

# ==========================================
# 4. EXPORTAR ORBITALES (Lo hacemos ANTES del ZFS para asegurar)
# ==========================================
print("\n Generando archivos .cube para VESTA...")
mo_alfa = mf.mo_coeff[0]
N_homo = mol.nelec[0] - 1

# Exportamos el HOMO (el más afectado por el campo)
cubegen.orbital(mol, 'NV_PASIVADO_HOMO.cube', mo_alfa[:, N_homo])
# Exportamos el HOMO-1
cubegen.orbital(mol, 'NV_PASIVADO_HOMO_1.cube', mo_alfa[:, N_homo - 1])
print("¡Archivos .cube creados con éxito!")

# ==========================================
# 5. ZFS (Solo para ver los números de nuevo)
# ==========================================
print("\n Calculando ZFS...")
calc_zfs = zfs.uhf.ZFS(mf)
try:
    calc_zfs.kernel()
except Exception:
    print("\n(Ignorando error final del logger. Revisa los valores D y E arriba)")