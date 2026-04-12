import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, steadystate, expect

# ==========================================================
# 1. DEFINICIÓN DE LA BASE Y OPERADORES
# ==========================================================
g0 = basis(7, 0); gm = basis(7, 1); gp = basis(7, 2)
e0 = basis(7, 3); em = basis(7, 4); ep = basis(7, 5)
s  = basis(7, 6)

# Proyector de fotoluminiscencia
P_excited = e0*e0.dag() + em*em.dag() + ep*ep.dag()

# Operadores de Espín Z (FALTABA Sz_es)
Sz_gs = 0*g0*g0.dag() - 1*gm*gm.dag() + 1*gp*gp.dag()
Sz_es = 0*e0*e0.dag() - 1*em*em.dag() + 1*ep*ep.dag()

# Operador de Espín X para el estado fundamental (acoplamiento de microondas)
Sx_gs = (1/np.sqrt(2)) * (g0*gm.dag() + gm*g0.dag() + g0*gp.dag() + gp*g0.dag())

# ==========================================================
# 2. PARÁMETROS FÍSICOS
# ==========================================================
Dg = 2.87e9         # ZFS GS (Hz)
De = 1.4e9          # ZFS ES (Hz) - FALTABA
gamma_e = 2.8024e6  # Ratio giromagnético (Hz/Gauss)

Bz = 50.0           # Campo estático a detectar (Gauss)
Omega_mw = 2e6      # Potencia de microondas (Frecuencia de Rabi en Hz)

# Tasas cinéticas de Lindblad (Hz)
kpump = 20e6; Orad = 80e6; kepms = 60e6; ke0s = 1e6; kgs0 = 3e6; kgspm = 0.2e6

c_ops = [
    np.sqrt(kpump)*e0*g0.dag(), np.sqrt(kpump)*em*gm.dag(), np.sqrt(kpump)*ep*gp.dag(),
    np.sqrt(Orad)*g0*e0.dag(),  np.sqrt(Orad)*gm*em.dag(),  np.sqrt(Orad)*gp*ep.dag(),
    np.sqrt(ke0s)*s*e0.dag(),   np.sqrt(kepms)*s*em.dag(),  np.sqrt(kepms)*s*ep.dag(),
    np.sqrt(kgs0)*g0*s.dag(),   np.sqrt(kgspm)*gm*s.dag(),  np.sqrt(kgspm)*gp*s.dag()
]

# ==========================================================
# 3. SIMULACIÓN CW-ODMR (BARRIDO DE FRECUENCIA)
# ==========================================================
# Barremos +/- 200 MHz alrededor del centro Dg
f_mw_list = np.linspace(Dg - 200e6, Dg + 200e6, 300)
pl_spectrum = []

print(f"Simulando espectro CW-ODMR para Bz = {Bz} Gauss...")

for f_mw in f_mw_list:
    # Estado Fundamental en el Marco Rotatorio (Restamos f_mw al ZFS)
    H_gs = (Dg - f_mw) * (Sz_gs**2) + gamma_e * Bz * Sz_gs + Omega_mw * Sx_gs
    
    # Estado Excitado en el marco del laboratorio (No afectado por la MW)
    H_es = De * (Sz_es**2) + gamma_e * Bz * Sz_es
    
    # Hamiltoniano Total (En frecuencia angular)
    H_total = 2 * np.pi * (H_gs + H_es) 
    
    # Resolver estado estacionario
    rho_ss = steadystate(H_total, c_ops) #Usar steadystate equivale a tener el laser aplicado continuamente
    pl_spectrum.append(expect(P_excited, rho_ss))

print("¡Simulación completada!")

# ==========================================================
# 4. GRAFICAR RESULTADOS
# ==========================================================
plt.figure(figsize=(9, 5))

# Normalizamos para que la PL fuera de resonancia sea ~1.0
pl_norm = np.array(pl_spectrum) / max(pl_spectrum)

plt.plot(f_mw_list / 1e9, pl_norm, color='navy', linewidth=2)
plt.xlabel('MW (GHz)', fontsize=12)
plt.ylabel('PL (Normalized)', fontsize=12)

# Líneas de referencia teóricas para los dos picos
plt.axvline(x=(Dg - gamma_e * Bz)/1e9, color='red', linestyle='--', alpha=0.5)
plt.axvline(x=(Dg + gamma_e * Bz)/1e9, color='red', linestyle='--', alpha=0.5)

plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()