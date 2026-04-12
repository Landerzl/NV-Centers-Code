import numpy as np
from pyscf import gto, scf
from pyscf.prop import zfs
from pyscf.scf import addons
from pyscf.tools import cubegen

# ==========================================
# 1. SISTEMA (1 NV Pasivado - Optimizado para RAM)
# ==========================================
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

mol = gto.M(
    atom = nv_pasivado,
    # AQUÍ ESTÁ LA MAGIA DE LA RAM: Base mixta
    basis = {
        'N': 'def2-svp',  # Alta resolución (Física importante)
        'C': 'def2-svp',  # Alta resolución (Física importante)
        'H': 'sto-3g'     # Baja resolución (Ahorro masivo de memoria)
    },
    charge = -1,
    spin = 2,
    verbose = 4
)

# Límite de RAM para PySCF (8000 MB = 8 GB). Ajusta según tu PC.
mol.max_memory = 8000 

# ==========================================
# 2. CAMPO ELÉCTRICO
# ==========================================
E_x, E_y, E_z = 0.002, 0.0, 0.0  

# ==========================================
# 3. CÁLCULO DEL ESTADO FUNDAMENTAL
# ==========================================
print("\n--- 1. Calculando Estado Fundamental ---")
mf_gs = scf.UKS(mol)
mf_gs.xc = 'b3lyp'

hcore_orig = mf_gs.get_hcore
def get_hcore_con_campo(mol=mol):
    h = hcore_orig(mol)
    mol.set_common_orig([0, 0, 0])
    dipole_ints = mol.intor_symmetric('int1e_r', comp=3)
    return h + (E_x * dipole_ints[0] + E_y * dipole_ints[1] + E_z * dipole_ints[2])

mf_gs.get_hcore = get_hcore_con_campo
mf_gs.kernel()

# ==========================================
# 4. EXCITACIÓN MANUAL (Delta-SCF)
# ==========================================
print("\n--- 2. Promoviendo electrón al Estado Excitado ---")
mo_occ_exc = mf_gs.mo_occ.copy()
mo_coeff = mf_gs.mo_coeff.copy()

N_beta = mol.nelec[1]
idx_a1 = N_beta - 1  # Índice del HOMO beta (Piso ocupado)
idx_e  = N_beta      # Índice del LUMO beta (Piso vacío)

mo_occ_exc[1][idx_a1] = 0.0  # Vaciamos el nivel de abajo
mo_occ_exc[1][idx_e]  = 1.0  # Forzamos el electrón arriba

# ==========================================
# 5. CÁLCULO DEL ESTADO EXCITADO (MOM)
# ==========================================
print("\n--- 3. Relajando el Estado Excitado (MOM) ---")
mf_exc = scf.UKS(mol)
mf_exc.xc = 'b3lyp'

# --- NUEVOS PARÁMETROS DE CONVERGENCIA ---
mf_exc.level_shift = 0.4  # Hacemos el "freno" más fuerte
mf_exc.max_cycle = 150    # Le damos 150 intentos en lugar de 50
# -----------------------------------------

mf_exc.get_hcore = get_hcore_con_campo

# El candado MOM para que no decaiga
mf_exc = addons.mom_occ(mf_exc, mo_coeff, mo_occ_exc)

dm_exc = mf_exc.make_rdm1(mo_coeff, mo_occ_exc)
mf_exc.kernel(dm0=dm_exc)
# ==========================================
# 6. EXPORTAR Y ZFS DEL ESTADO EXCITADO
# ==========================================
print("\n--- 4. Resultados del Estado Excitado ---")
cubegen.orbital(mol, 'EXCITADO_hueco_a1.cube', mf_exc.mo_coeff[1][:, idx_a1])
cubegen.orbital(mol, 'EXCITADO_electron_e.cube', mf_exc.mo_coeff[1][:, idx_e])
print("-> Archivos .cube generados. Ábrelos en VESTA.")

print("\n Calculando ZFS (Efecto Stark Excitado)...")
calc_zfs = zfs.uhf.ZFS(mf_exc)
calc_zfs.verbose = 4  

try:
    calc_zfs.kernel()
except TypeError:
    pass # Ignoramos el bug de impresión del logger
except Exception as e:
    print(f"\n[AVISO] Hubo un error interno en el ZFS: {e}")

print("\n" + "="*45)
print("   RESULTADOS FINALES (ESTADO EXCITADO)   ")
print("="*45)
try:
    print(f"Parámetro Axial (D):     {calc_zfs.D_value:.6f} cm^-1")
    print(f"Parámetro Rómbico (E):   {calc_zfs.E_value:.6f} cm^-1")
    print("\n*Fíjate si el parámetro E es mayor que en el estado fundamental.*")
except AttributeError:
    print("❌ No se pudieron extraer D y E.")
    print("Si el SCF no convergió, sube el valor de 'level_shift' a 0.3 o 0.4.")
print("="*45)