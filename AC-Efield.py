import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, steadystate, mesolve, expect

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
gamma_e = 2.8024e6  # Hz/Gauss

dgspara = 0.35      # Hz / (V/cm)
despara = 20e3      # Hz / (V/cm)

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
# 3. CONFIGURACIÓN DEL EXPERIMENTO DINÁMICO
# ==========================================================
# Fijamos un ángulo pequeño (1.5 grados) para un valle ESLAC afilado
alpha = np.radians(1.5) 

# B_bias fijo en la zona de máxima pendiente del ESLAC (~502 Gauss)
B_bias = 502.0
Bz = B_bias * np.cos(alpha)
Bx = B_bias * np.sin(alpha)

# Frecuencia del campo eléctrico oscilante (1 MHz)
f_osc = 1e6
omega = 2 * np.pi * f_osc

# PARTE ESTÁTICA DEL HAMILTONIANO (H0)
H0_gs = Dg * (Sz_gs**2) + gamma_e * Bz * Sz_gs + gamma_e * Bx * Sx_gs
H0_es = De * (Sz_es**2) + gamma_e * Bz * Sz_es + gamma_e * Bx * Sx_es
H0 = 2 * np.pi * (H0_gs + H0_es)

# PARTE DEPENDIENTE DEL CAMPO ELÉCTRICO (H1) - Efecto Stark
# El operador H1 se multiplicará por la función E(t)
H1_stark_gs = dgspara * (Sz_gs**2)
H1_stark_es = despara * (Sz_es**2)
H1 = 2 * np.pi * (H1_stark_gs + H1_stark_es)

# Función de coeficiente dependiente del tiempo: E(t) = E0 * sin(omega * t)
def H1_coeff(t, args):
    return args['E0'] * np.sin(args['omega'] * t)

# ==========================================================
# 4. SIMULACIÓN CON MESOLVE (Evolución Temporal)
# ==========================================================
# Simularemos 3 periodos completos del campo eléctrico
tiempo_total = 3 / f_osc
tlist = np.linspace(0, tiempo_total, 400)

# Estado inicial (rho0): Encontramos el equilibrio estático en t=0 (sin E-field)
# Esto evita ver transitorios de encendido del láser, solo vemos la respuesta al campo
rho0 = steadystate(H0, c_ops)

def simular_senal(E0_val):
    print(f"Simulando para E0 = {E0_val} V/cm...")
    args = {'E0': E0_val, 'omega': omega}
    
    # Formato QuTiP para Hamiltoniano dependiente del tiempo
    H_t = [H0, [H1, H1_coeff]]
    
    # mesolve integra la ecuación de Lindblad en el tiempo
    resultado = mesolve(H_t, rho0, tlist, c_ops, [P_excited], args=args)
    return resultado.expect[0]

# Simulamos para un campo débil y un campo fuerte
pl_lineal = simular_senal(200)
pl_saturado = simular_senal(1500)

# ==========================================================
# 5. GRAFICAR LOS RESULTADOS (OSCILOSCOPIO VIRTUAL)
# ==========================================================
# Creamos la señal de referencia del campo eléctrico (normalizada) para comparar la fase
referencia_E = np.sin(omega * tlist)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Panel 1: Régimen Lineal
ax1.plot(tlist * 1e6, pl_lineal, color='blue', linewidth=2, label='Señal PL (NV)')
ax1.plot(tlist * 1e6, referencia_E * (np.ptp(pl_lineal)/2) + np.mean(pl_lineal), 
         color='grey', linestyle='--', alpha=0.5, label='Fase del Campo $E(t)$')
ax1.set_title('($E_0 = 200$ V/cm)', fontsize=12)
ax1.set_ylabel('PL (u.a.)')
ax1.grid(True, alpha=0.3)
ax1.legend(loc='upper right')

# Panel 2: Régimen Saturado / Distorsionado
ax2.plot(tlist * 1e6, pl_saturado, color='red', linewidth=2, label='Señal PL (NV)')
ax2.plot(tlist * 1e6, referencia_E * (np.ptp(pl_saturado)/2) + np.mean(pl_saturado), 
         color='grey', linestyle='--', alpha=0.5, label='Fase del Campo $E(t)$')
ax2.set_title('($E_0 = 1500$ V/cm)', fontsize=12)
ax2.set_xlabel('Tiempo ($\mu s$)', fontsize=12)
ax2.set_ylabel('PL (u.a.)')
ax2.grid(True, alpha=0.3)
ax2.legend(loc='upper right')

plt.tight_layout()
plt.show()