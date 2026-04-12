import numpy as np
from pyscf import gto, scf
from pyscf.prop import zfs
from pyscf.scf import addons
import matplotlib.pyplot as plt
import csv

# Importaciones para el hackeo de la terminal
import io
import re
from contextlib import redirect_stdout

# ==========================================
# 1. SISTEMA (1 NV Pasivado - Base Mixta)
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
    basis = {'N': 'def2-svp', 'C': 'def2-svp', 'H': 'sto-3g'},
    charge = -1,
    spin = 2,
    verbose = 3
)
mol.max_memory = 8000 

# ==========================================
# 2. DEFINIR EL BARRIDO DE CAMPO ELÉCTRICO
# ==========================================
# Hacemos 5 puntos desde 0 hasta 0.002 a.u.
E_x_vals = np.linspace(0.0, 0.002, 5)
E_y, E_z = 0.0, 0.0

campos_exitosos = []
valores_D = []
valores_E = []

print(f"\nINICIANDO BARRIDO STARK EXCITADO: {len(E_x_vals)} puntos.")
print("="*50)

# ==========================================
# 3. BUCLE DE CÁLCULO
# ==========================================
for E_x in E_x_vals:
    campo_kv_cm = E_x * (5.14e6) 
    print(f"\n>>> Calculando para E_x = {E_x:.5f} a.u. (~{campo_kv_cm:.0f} kV/cm)")
    
    try: # <--- ESTE ES EL TRY QUE SE QUEDÓ HUÉRFANO ANTES
        # --- A. FUNDAMENTAL ---
        mf_gs = scf.UKS(mol)
        mf_gs.xc = 'b3lyp'
        
        hcore_orig = mf_gs.get_hcore
        def get_hcore_con_campo(mol=mol, ex=E_x):
            h = hcore_orig(mol)
            mol.set_common_orig([0, 0, 0])
            dipole_ints = mol.intor_symmetric('int1e_r', comp=3)
            return h + (ex * dipole_ints[0] + E_y * dipole_ints[1] + E_z * dipole_ints[2])
            
        mf_gs.get_hcore = get_hcore_con_campo
        mf_gs.kernel()
        
        # --- B. EXCITACIÓN ---
        mo_occ_exc = mf_gs.mo_occ.copy()
        mo_coeff = mf_gs.mo_coeff.copy()
        N_beta = mol.nelec[1]
        idx_a1 = N_beta - 1
        idx_e  = N_beta
        mo_occ_exc[1][idx_a1] = 0.0
        mo_occ_exc[1][idx_e]  = 1.0
        
        # --- C. EXCITADO CON MOM ---
        mf_exc = scf.UKS(mol)
        mf_exc.xc = 'b3lyp'
        mf_exc.level_shift = 0.4
        mf_exc.max_cycle = 150
        mf_exc.get_hcore = get_hcore_con_campo
        
        mf_exc = addons.mom_occ(mf_exc, mo_coeff, mo_occ_exc)
        dm_exc = mf_exc.make_rdm1(mo_coeff, mo_occ_exc)
        mf_exc.kernel(dm0=dm_exc)
        
        if not mf_exc.converged:
            print(f" [!] Aviso: SCF excitado inestable. Forzando extracción ZFS...")
            
        # --- D. ZFS (EXTRACCIÓN EXTREMA) ---
        mf_exc.converged = True # El engaño principal
        
        calc_zfs = zfs.uhf.ZFS(mf_exc)
        calc_zfs.verbose = 4
        
        # La trampa de la terminal
        trampa = io.StringIO()
        with redirect_stdout(trampa):
            try:
                calc_zfs.kernel()
            except Exception:
                pass # Ignoramos el crasheo interno del ZFS
                
        texto_terminal = trampa.getvalue()
        
        # El hackeo con expresiones regulares
        match_D = re.search(r"Axial\s+parameter D =\s*([+-]?\d+\.\d+)", texto_terminal)
        match_E = re.search(r"Rhombic\s+parameter E =\s*([+-]?\d+\.\d+)", texto_terminal)
        
        if match_D and match_E:
            D_val = float(match_D.group(1))
            E_val = float(match_E.group(1))
            
            campos_exitosos.append(campo_kv_cm)
            valores_D.append(D_val)
            valores_E.append(E_val)
            
            print(f" [+] Datos capturados a la fuerza -> D: {D_val:.5f} | E: {E_val:.5f}")
        else:
            print(" [X] PySCF colapsó y no se imprimieron los números.")

    # ESTE ES EL EXCEPT QUE FALTABA
    except Exception as error:
        print(f" [X] Error inesperado en el bucle principal (E_x={E_x}): {error}")

print("\n" + "="*50)
print("BARRIDO COMPLETADO.")

# ==========================================
# 4. GUARDAR DATOS Y GRAFICAR
# ==========================================
if len(campos_exitosos) > 0:
    with open('resultados_zfs_excitado.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Campo (kV/cm)', 'D (cm^-1)', 'E (cm^-1)'])
        for i in range(len(campos_exitosos)):
            writer.writerow([campos_exitosos[i], valores_D[i], valores_E[i]])
    print("Datos guardados en 'resultados_zfs_excitado.csv'.")

    plt.figure(figsize=(10, 5))

    plt.subplot(1, 2, 1)
    plt.plot(campos_exitosos, valores_D, marker='o', color='blue', linestyle='dashed')
    plt.title('Parámetro D (Estado Excitado)')
    plt.xlabel('Campo Eléctrico $E_x$ (kV/cm)')
    plt.ylabel('ZFS Axial D (cm$^{-1}$)')
    plt.grid(True)

    plt.subplot(1, 2, 2)
    plt.plot(campos_exitosos, valores_E, marker='o', color='red', linestyle='dashed')
    plt.title('Parámetro E (Efecto Stark Rómbico)')
    plt.xlabel('Campo Eléctrico $E_x$ (kV/cm)')
    plt.ylabel('ZFS Rómbico E (cm$^{-1}$)')
    plt.grid(True)

    plt.tight_layout()
    plt.savefig('Efecto_Stark_ZFS_Excitado.png', dpi=300)
    print("¡Gráfica generada! Abre 'Efecto_Stark_ZFS_Excitado.png'.")
else:
    print("No se pudo extraer ningún dato.")