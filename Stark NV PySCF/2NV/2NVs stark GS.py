import numpy as np
from pyscf import gto, scf
from pyscf.prop import zfs
from pyscf.tools import cubegen

# ==========================================
# 1. SISTEMA (Dos centros NV enfrentados)
# ==========================================
# Hemos diseñado un clúster con dos Nitrógenos (N1 y N2) 
# y sus correspondientes vacantes.
dos_nvs = '''
N    0.000000    0.000000    0.000000
C    1.020000    1.020000    1.020000
C   -1.020000   -1.020000    1.020000
C    1.020000   -1.020000   -1.020000
N    0.000000    0.000000    3.500000
C    1.020000    1.020000    2.480000
C   -1.020000   -1.020000    2.480000
C    1.020000   -1.020000    4.520000
H    1.500000    2.000000    0.800000
H   -1.500000   -2.000000    0.800000
H    1.500000    2.000000    2.700000
H   -1.500000   -2.000000    2.700000
'''

mol = gto.M(
    atom = dos_nvs,
    basis = 'sto-3g', #esta base es una mierda, pero me quedo sin RAM, el resultado es malísimo
    charge = -2,      # Dos cargas negativas (una por cada NV)
    spin = 4,         # Quintuplete (2 electrones de cada NV alineados = 4)
    verbose = 4
)

# ==========================================
# 2. CAMPO ELÉCTRICO Y CÁLCULO
# ==========================================
E_x, E_y, E_z = 0.002, 0.0, 0.0

mf = scf.UKS(mol)
mf.xc = 'b3lyp'
mf.level_shift = 0.3 # Más estabilidad para clúster grande

hcore_orig = mf.get_hcore
def get_hcore_con_campo(mol=mol):
    h = hcore_orig(mol)
    mol.set_common_orig([0, 0, 0])
    dipole_ints = mol.intor_symmetric('int1e_r', comp=3)
    return h + (E_x * dipole_ints[0] + E_y * dipole_ints[1] + E_z * dipole_ints[2])

mf.get_hcore = get_hcore_con_campo
mf.kernel()

# ==========================================
# 3. EXPORTAR Y ZFS (Versión Debug Extremo)
# ==========================================
print("\n[1/2] Generando orbitales...")
mo_alfa = mf.mo_coeff[0]
N_homo = mol.nelec[0] - 1
cubegen.orbital(mol, 'DOS_NV_HOMO.cube', mo_alfa[:, N_homo])
print(" -> Archivo DOS_NV_HOMO.cube listo.")

print("\n[2/2] Calculando ZFS (Modo Debug)...")
calc_zfs = zfs.uhf.ZFS(mf)

# CONFIGURACIÓN DE EMERGENCIA
calc_zfs.verbose = 9      # Máximo nivel de detalle en la terminal
calc_zfs.cphf = False     # Apagamos CPHF para evitar que el sistema lineal colapse
# Si hiciste el parche en uhf.py (zfs_soc = 0.0), esto debería fluir.

# Ejecutamos SIN el bloque try-except para ver el error real si es que truena
zfs_tensor = calc_zfs.kernel()

# Si llega aquí, imprimimos manual
print(f"\n--- VALORES EXTRAÍDOS ---")
print(f"D = {calc_zfs.D_value} cm^-1")
print(f"E = {calc_zfs.E_value} cm^-1")