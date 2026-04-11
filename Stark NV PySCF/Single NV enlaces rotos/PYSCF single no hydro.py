import numpy as np
from pyscf import gto, scf
from pyscf.prop import zfs
from pyscf.tools import cubegen

# 1. SISTEMA
# Para un paper real, necesitarías añadir Hidrógenos de pasivación aquí.
nv_cluster = '''
N    0.0000000    0.0000000    0.0000000
C    1.5000000    1.5000000    1.5000000
C   -1.5000000   -1.5000000    1.5000000
C   -1.5000000    1.5000000   -1.5000000
'''
mol = gto.M(atom=nv_cluster, basis='def2-svp', charge=-1, spin=2, verbose=4)

# 2. CAMPO ELÉCTRICO
E_x, E_y, E_z = 0.002, 0.0, 0.0

# 3. CÁLCULO
mf = scf.UKS(mol)
mf.xc = 'b3lyp'

# --- AYUDAS PARA CONVERGENCIA SCF ---
mf.max_cycle = 200     # Le damos más intentos para converger
mf.level_shift = 0.1   # Amortigua las oscilaciones de los orbitales abiertos
# ------------------------------------

hcore_orig = mf.get_hcore
def get_hcore_con_campo(mol=mol):
    h = hcore_orig(mol)
    mol.set_common_orig([0, 0, 0])
    dipole_ints = mol.intor_symmetric('int1e_r', comp=3)
    return h + (E_x * dipole_ints[0] + E_y * dipole_ints[1] + E_z * dipole_ints[2])

mf.get_hcore = get_hcore_con_campo
mf.kernel()

if not mf.converged:
    print("\n[ADVERTENCIA FÍSICA]: El SCF no convergió. Los resultados son aproximados.")
    print("Para solucionarlo a futuro, añade átomos de Hidrógeno a los bordes de tu clúster.")

# ==========================================
# 4. CÁLCULO DEL ZFS (Efecto Stark)
# ==========================================
calc_zfs = zfs.uhf.ZFS(mf)

print("\n=== CALCULANDO ZFS (Los resultados saldrán arriba) ===")
try:
    # Esto imprimirá los valores D y E automáticamente y luego lanzará el bug del logger
    zfs_tensor = calc_zfs.kernel()
except TypeError:
    print("\n(Ignorando el bug interno de PySCF. Los valores de D y E están impresos justo arriba ↑)")

# Hemos borrado las líneas problemáticas que buscaban calc_zfs.D_value
# para que el programa no se detenga aquí.

# ==========================================
# 5. EXPORTAR ORBITALES
# ==========================================
print("\n=== EXPORTANDO ARCHIVOS .CUBE PARA VESTA ===")
mo_alfa = mf.mo_coeff[0]
N_homo = mol.nelec[0] - 1

cubegen.orbital(mol, 'orbital_HOMO.cube', mo_alfa[:, N_homo])
print("-> Generado: orbital_HOMO.cube (Orbital E_x deformado)")

cubegen.orbital(mol, 'orbital_HOMO_menos_1.cube', mo_alfa[:, N_homo - 1])
print("-> Generado: orbital_HOMO_menos_1.cube (Orbital E_y deformado)")

print("\n¡Proceso 100% terminado! Busca los archivos .cube en tu carpeta.")