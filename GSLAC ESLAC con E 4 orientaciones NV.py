import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, steadystate, expect

# ==========================================================
# 1. DEFINICIÓN DE LA BASE Y OPERADORES
# ==========================================================
g0 = basis(7, 0); gm = basis(7, 1); gp = basis(7, 2)
e0 = basis(7, 3); em = basis(7, 4); ep = basis(7, 5)
s  = basis(7, 6)

P_excited = e0*e0.dag() + em*em.dag() + ep*ep.dag()

Sz_gs = 0*g0*g0.dag() - 1*gm*gm.dag() + 1*gp*gp.dag()
Sz_es = 0*e0*e0.dag() - 1*em*em.dag() + 1*ep*ep.dag()

Sx_gs = (1/np.sqrt(2)) * (g0*gm.dag() + gm*g0.dag() + g0*gp.dag() + gp*g0.dag())
Sx_es = (1/np.sqrt(2)) * (e0*em.dag() + em*e0.dag() + e0*ep.dag() + ep*e0.dag())

# ==========================================================
# 2. PARÁMETROS DEL SISTEMA
# ==========================================================
Dg = 2.87e9       # ZFS GS (Hz)
De = 1.4e9        # ZFS ES (Hz)
gamma_e = 2.8024e6  # Ratio giromagnético (Hz/Gauss)

# Parámetros Stark longitudinales (Hz / (V/cm))
dgspara = 0.35
despara = 20e3

# Tasas cinéticas (Hz)
kpump = 20e6; Orad = 80e6
kepms = 60e6; ke0s = 1e6
kgs0  = 3e6;  kgspm = 0.2e6

c_ops = [
    np.sqrt(kpump)*e0*g0.dag(), np.sqrt(kpump)*em*gm.dag(), np.sqrt(kpump)*ep*gp.dag(),
    np.sqrt(Orad)*g0*e0.dag(),  np.sqrt(Orad)*gm*em.dag(),  np.sqrt(Orad)*gp*ep.dag(),
    np.sqrt(ke0s)*s*e0.dag(),   np.sqrt(kepms)*s*em.dag(),  np.sqrt(kepms)*s*ep.dag(),
    np.sqrt(kgs0)*g0*s.dag(),   np.sqrt(kgspm)*gm*s.dag(),  np.sqrt(kgspm)*gp*s.dag()
]

# ==========================================================
# 3. GEOMETRÍA DEL CRISTAL Y LOS CAMPOS EXTERNOS
# ==========================================================
n1 = np.array([ 1,  1,  1]) / np.sqrt(3)
n2 = np.array([-1, -1,  1]) / np.sqrt(3)
n3 = np.array([ 1, -1, -1]) / np.sqrt(3)
n4 = np.array([-1,  1, -1]) / np.sqrt(3)
orientaciones_NV = [n1, n2, n3, n4]

# Ángulos esféricos EXACTOS para el eje n1 [111]
theta_n1 = np.arccos(1/np.sqrt(3))  # ~54.74 grados
phi_n1 = np.radians(45)             # 45 grados

# Añadimos una pequeñísima desalineación (ej. 1.5 grados) 
# Esto genera el B_x justo para ver el valle sin destruirlo.
desalineacion = np.radians(1.5) 

theta_B = theta_n1 + desalineacion
phi_B   = phi_n1

# Vector de campo magnético casi alineado con n1
u_B = np.array([
    np.sin(theta_B) * np.cos(phi_B),
    np.sin(theta_B) * np.sin(phi_B),
    np.cos(theta_B)
])

# Dirección del Campo Eléctrico aplicado [001]
u_E = np.array([0, 0, 1])

# ==========================================================
# 4. FUNCIÓN DEL HAMILTONIANO (Dependiente de orientaciones)
# ==========================================================
def build_hamiltonian_nv(B_mag, E_mag_lab, nv_dir):
    # Vectores totales en el marco del laboratorio
    B_vec = B_mag * u_B
    E_vec = E_mag_lab * u_E
    
    # Proyecciones sobre el eje de este NV específico
    Bz_eff = np.dot(B_vec, nv_dir)
    Ez_eff = np.dot(E_vec, nv_dir)
    
    # Campo magnético transversal efectivo (mezcla cuántica)
    B_trans_vec = np.cross(B_vec, nv_dir)
    Bx_eff = np.linalg.norm(B_trans_vec)
    
    # Efecto Stark: El ZFS cambia según la proyección del campo eléctrico
    D_gs_eff = Dg + dgspara * Ez_eff
    D_es_eff = De + despara * Ez_eff
    
    # Construcción del Hamiltoniano
    H_gs_diag = D_gs_eff * (Sz_gs**2) + gamma_e * Bz_eff * Sz_gs
    H_es_diag = D_es_eff * (Sz_es**2) + gamma_e * Bz_eff * Sz_es
    
    H_gs_mix = gamma_e * Bx_eff * Sx_gs
    H_es_mix = gamma_e * Bx_eff * Sx_es
    
    H_Hz = (H_gs_diag + H_gs_mix) + (H_es_diag + H_es_mix)
    return 2 * np.pi * H_Hz

# ==========================================================
# 5. SIMULACIÓN (BARRIDO B Y E)
# ==========================================================
# Rango B para capturar los ESLACs desplazados y los GSLACs
B_vals = np.linspace(400, 1600, 300) 

# Campos eléctricos macroscópicos a simular (V/cm)
Ez_lab_values = [5e3, 10e3, 30e3]
resultados_pl = {}

print("Iniciando simulaciones masivas (Ensamble + Stark)...")

for E_lab in Ez_lab_values:
    print(f"-> Calculando para E_lab = {E_lab/1000} kV/cm...")
    pl_ensamble = []
    
    for B in B_vals:
        pl_total = 0
        # Promediamos sobre las 4 orientaciones para un B y E dados
        for nv_dir in orientaciones_NV:
            H = build_hamiltonian_nv(B, E_lab, nv_dir)
            rho_ss = steadystate(H, c_ops)
            pl_total += expect(P_excited, rho_ss)
        
        pl_ensamble.append(pl_total / 4)
        
    # Guardamos y normalizamos
    resultados_pl[E_lab] = np.array(pl_ensamble) / max(pl_ensamble)

print("¡Simulaciones completadas!")

# ==========================================================
# 6. GRAFICAR RESULTADOS
# ==========================================================
plt.figure(figsize=(12, 7))
colores = ['#1f77b4', '#ff7f0e', '#d62728'] # Azul, Naranja, Rojo

for i, E_lab in enumerate(Ez_lab_values):
    etiqueta = f'$E_{{lab}}$ = {int(E_lab/1000)} kV/cm'
    # Usamos un offset (desplazamiento vertical) para separar las curvas y verlas mejor
    offset = -i * 0.15 
    plt.plot(B_vals, resultados_pl[E_lab] + offset, color=colores[i], linewidth=2, label=etiqueta)


plt.xlabel('$|B|$ (Gauss)', fontsize=12)
plt.ylabel('PL (Apilada)', fontsize=12)
plt.yticks([]) # Ocultamos los números del eje Y porque están apiladas
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()