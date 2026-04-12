import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, mesolve, expect, Options

# ==========================================================
# 1. DEFINICIÓN DE LA BASE Y OPERADORES (Subespacio GS)
# ==========================================================
# Aunque es un sistema de 7 niveles, para la dinámica coherente de Ramsey
# a oscuras, solo necesitamos operar sobre el estado fundamental.
g0 = basis(7, 0)
gm = basis(7, 1)
gp = basis(7, 2)

# Operador de medición (Población en el estado brillante m_s = 0)
P_0 = g0 * g0.dag()

# Operadores de Espín Z y X para el estado fundamental
Sz_gs = 0*g0*g0.dag() - 1*gm*gm.dag() + 1*gp*gp.dag()
Sx_gs = (1/np.sqrt(2)) * (g0*gm.dag() + gm*g0.dag() + g0*gp.dag() + gp*g0.dag())

# ==========================================================
# 2. PARÁMETROS DEL EXPERIMENTO
# ==========================================================
Dg = 2.87e9         # ZFS (Hz)
gamma_e = 2.8024e6  # Ratio giromagnético (Hz/Gauss)
Bz = 50.0           # Campo estático que queremos medir (Gauss)

# Frecuencia de resonancia exacta
f_res = Dg - gamma_e * Bz

# Configuración del Pulso de Microondas
delta = 2e6         # Desafinación de 2 MHz (para generar las franjas)
f_mw = f_res + delta
Omega_mw = 10e6     # Potencia MW (10 MHz)

# Cálculo del tiempo del pulso pi/2
# (Para espín 1, la frecuencia efectiva de Rabi es sqrt(2)*Omega)
Omega_eff = np.sqrt(2) * Omega_mw
t_pi2 = 1 / (4 * Omega_eff) 

# Coherencia T2* (Microsegundos)
T2_star = 1.0e-6    # 1 us de tiempo de coherencia transversal

# Operador de colapso de PURO DESFASE (Pérdida de coherencia)
c_ops_ramsey = [np.sqrt(1 / T2_star) * Sz_gs]

# ==========================================================
# 3. CONSTRUCCIÓN DE LOS HAMILTONIANOS (Marco Rotatorio)
# ==========================================================
# Hamiltoniano a oscuras (Solo Zeeman y ZFS efectivo)
H_free = 2 * np.pi * ((Dg - f_mw) * (Sz_gs**2) + gamma_e * Bz * Sz_gs)

# Hamiltoniano durante el pulso (Se añade la interacción transversal con MW)
H_pulse = H_free + 2 * np.pi * (Omega_mw * Sx_gs)

# ==========================================================
# 4. SIMULACIÓN DE LA SECUENCIA DE RAMSEY
# ==========================================================
# Vector de tiempos de evolución libre (hasta 3 microsegundos)
tau_list = np.linspace(0, 3e-6, 150)
ramsey_signal = []

# Estado inicial: Inicialización perfecta por láser en m_s = 0
rho_init = g0 * g0.dag()

# ¡LA CORRECCIÓN CLAVE! Aumentamos los pasos permitidos para el integrador
opciones = Options(nsteps=50000)

print(f"Simulando Secuencia Ramsey (Desafinación = {delta/1e6} MHz)...")

for tau in tau_list:
    # --- Bloque 1: Primer Pulso pi/2 ---
    res_pulse1 = mesolve(H_pulse, rho_init, [0, t_pi2], c_ops_ramsey, options=opciones)
    rho_1 = res_pulse1.states[-1] 
    
    # --- Bloque 2: Evolución Libre (Tiempo tau) a oscuras ---
    if tau == 0:
        rho_2 = rho_1
    else:
        res_free = mesolve(H_free, rho_1, [0, tau], c_ops_ramsey, options=opciones)
        rho_2 = res_free.states[-1]
        
    # --- Bloque 3: Segundo Pulso pi/2 ---
    res_pulse2 = mesolve(H_pulse, rho_2, [0, t_pi2], c_ops_ramsey, options=opciones)
    rho_final = res_pulse2.states[-1]
    
    # Medida final: Extraemos la probabilidad de encontrar el espín en |0>
    ramsey_signal.append(expect(P_0, rho_final))

print("¡Simulación completada con éxito!")

# ==========================================================
# 5. GRAFICAR RESULTADOS
# ==========================================================
plt.figure(figsize=(9, 5))

# Convertimos el tiempo a microsegundos para la gráfica
plt.plot(tau_list * 1e6, ramsey_signal, color='purple', linewidth=2)

plt.title(f'Interferometría de Ramsey ($\delta$ = {delta/1e6} MHz, $T_2^*$ = {T2_star*1e6} $\mu$s)', fontsize=14)
plt.xlabel('Tiempo de Evolución Libre $\tau$ ($\mu$s)', fontsize=12)
plt.ylabel('Población en estado $|0\\rangle$', fontsize=12)

# Ajustamos los límites Y para que la oscilación se vea centrada entre 0 y 1
plt.ylim(0, 1)

plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()