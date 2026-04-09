import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, steadystate, expect

# ==========================================================
# 1. DEFINICIÓN DE LA BASE Y OPERADORES DE ESPÍN
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
Dg = 2.87e9
De = 1.4e9
gamma_e = 2.8024e6  # Hz / Gauss

# Parámetros Stark longitudinales (Hz / (V/cm))
dgspara = 0.35
despara = 20e3

# Ángulo de desalineación del campo B
alpha = 0.1 * (np.pi / 2)

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
# 3. FUNCIÓN DEL HAMILTONIANO (Dependiente de B y Ez)
# ==========================================================
def build_hamiltonian(B_mag, Ez):
    Bz = B_mag * np.cos(alpha)
    Bx = B_mag * np.sin(alpha)
    
    # Efecto Stark: Modificación del ZFS por el campo eléctrico axial
    D_gs_eff = Dg + dgspara * Ez
    D_es_eff = De + despara * Ez
    
    H_gs_diag = D_gs_eff * (Sz_gs**2) + gamma_e * Bz * Sz_gs
    H_es_diag = D_es_eff * (Sz_es**2) + gamma_e * Bz * Sz_es
    
    H_gs_mix = gamma_e * Bx * Sx_gs
    H_es_mix = gamma_e * Bx * Sx_es
    
    H_Hz = (H_gs_diag + H_gs_mix) + (H_es_diag + H_es_mix)
    return 2 * np.pi * H_Hz

# ==========================================================
# 4. SIMULACIÓN PARA MÚLTIPLES CAMPOS ELÉCTRICOS
# ==========================================================
# Rango B ajustado para capturar el gran desplazamiento del ESLAC
B_vals = np.linspace(450, 1150, 300) 

# Valores de campo eléctrico a simular (en V/cm)
Ez_values = [5e3, 10e3, 30e3]

# Diccionario para guardar los resultados
resultados_pl = {}

print("Iniciando simulaciones...")
for Ez in Ez_values:
    print(f"Simulando para Ez = {Ez/1000} kV/cm...")
    pl_signal = []
    for B in B_vals:
        H = build_hamiltonian(B, Ez)
        rho_ss = steadystate(H, c_ops)
        pl_signal.append(expect(P_excited, rho_ss))
    
    # Normalizar cada señal por su propio máximo
    resultados_pl[Ez] = np.array(pl_signal) / max(pl_signal)

print("¡Simulaciones completadas!")

# ==========================================================
# 5. GRAFICAR RESULTADOS
# ==========================================================
plt.figure(figsize=(10, 6))

colores = ['#1f77b4', '#ff7f0e', '#d62728'] # Azul, Naranja, Rojo

for i, Ez in enumerate(Ez_values):
    etiqueta = f'$E_z$ = {int(Ez/1000)} kV/cm'
    plt.plot(B_vals, resultados_pl[Ez], color=colores[i], linewidth=2, label=etiqueta)

plt.axvline(x=1024, color='black', linestyle='--', alpha=0.5, label='GSLAC (Inmóvil)')

plt.title('Desplazamiento Stark del ESLAC vs Campo Eléctrico Axial', fontsize=14)
plt.xlabel('Magnitud del Campo Magnético B (Gauss)', fontsize=12)
plt.ylabel('Fotoluminiscencia Estacionaria (Norm.)', fontsize=12)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()