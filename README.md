# NV Center Simulations

Quantum open-system and Density Functional Theory (DFT) simulations of **Nitrogen-Vacancy (NV) centres in diamond**. This repository combines:

1. **Open-system dynamics** built with [QuTiP](https://qutip.org/), reproducing key spectroscopic observables (photoluminescence, Ramsey fringes, GSLAC/ESLAC features) as a function of applied magnetic and electric fields.
2. **Quantum chemistry calculations** using [PySCF](https://pyscf.org/) to evaluate the Stark effect, Zero-Field Splitting (ZFS) parameters (D and E), and calculate molecular orbitals representation (HOMO/LUMO) under different conditions (passivated NV, excited states, 2-NV interactions).

---

## Repository structure

```text
NV Centers Code/
├── 7 level Model + Lindbald/                    # QuTiP open-system dynamics
│   ├── AC-Efield.py                             # Dynamic response to an oscillating AC electric field
│   ├── CW-ODMR.py                               # Continuous-wave ODMR spectrum simulation
│   ├── GSLAC ESLAC con E.py                     # Static E-field shift of GSLAC and ESLAC
│   ├── GSLAC ESLAC con E 4 orientaciones NV.py  # 4-orientation ensemble average
│   ├── Ramsey.py                                # Pulsed Ramsey interferometry simulation
│   └── NV_E_G_SLAC_explanation.pdf              # Supplementary notes on GSLAC / ESLAC physics
├── Stark NV PySCF/                              # Quantum Chemistry and DFT calculations
│   ├── Single NV enlaces rotos/                 # Unpassivated single NV model
│   ├── Single NV pasivizado con H/              # H-passivated single NV model
│   ├── Excited state/                           # Delta-SCF and MOM computations for the excited state
│   └── 2NV pasivizado con H/                    # Model of two interacting NV centres
├── 2NV/                                         # Alternative / Experimental 2-NV code
└── Figures/                                     # Pre-generated output figures (PNG)
```

---

## Physics background

### 1. QuTiP 7-Level Model
The QuTiP scripts model the full **7-level system** (|g₀⟩, |g₋⟩, |g₊⟩, |e₀⟩, |e₋⟩, |e₊⟩, |s⟩) under continuous laser driving and incoherent decay (Lindblad formalism). The key parameters shared across all scripts are:

| Symbol | Value | Meaning |
|---|---|---|
| D_g | 2.87 GHz | Zero-field splitting, ground state |
| D_e | 1.4 GHz | Zero-field splitting, excited state |
| γ_e | 2.8024 MHz / Gauss | Electron gyromagnetic ratio |
| d_∥^GS | 0.35 Hz/(V/cm) | Axial Stark coefficient, GS |
| d_∥^ES | 20 kHz/(V/cm) | Axial Stark coefficient, ES |

### 2. PySCF DFT Calculations
Using the `pyscf.prop.zfs` module, we evaluate the D and E parameters directly from the spin density of the ground state (and excited state via Maximum Overlap Method, MOM). Additionally, `.cube` files are generated containing the spatial distribution of the molecular orbitals (e.g., HOMO) and can be visualized in VESTA.

---

## Requirements

```text
python >= 3.9
numpy
matplotlib
qutip >= 4.7
pyscf
```

Install with:

```bash
pip install numpy matplotlib qutip pyscf
```

> **Note:** QuTiP 5 introduced API changes. The scripts are tested against **QuTiP 4.7.x**. PySCF is used for the quantum chemistry calculations in the `Stark NV PySCF` folder.

---

## Running the simulations

### QuTiP dynamics (in `7 level Model + Lindbald`)
```bash
python "7 level Model + Lindbald/CW-ODMR.py"
python "7 level Model + Lindbald/Ramsey.py"
python "7 level Model + Lindbald/GSLAC ESLAC con E.py"
python "7 level Model + Lindbald/GSLAC ESLAC con E 4 orientaciones NV.py"
python "7 level Model + Lindbald/AC-Efield.py"
```

### PySCF calculations (in `Stark NV PySCF` subdirectories)
```bash
python "Stark NV PySCF/Single NV pasivizado con H/NV stark single + hydro.py"
python "Stark NV PySCF/Excited state/Estado excitado.py"
# ...and others
```
> **Performance Note:** PySCF calculations require significant RAM, especially for the 2NV model or high-resolution basis sets (`def2-svp`). You may need to tune `mol.max_memory` in the scripts to fit your hardware. Once `.cube` files are generated, use a software like **VESTA** to render the electron densities.

---

## Output gallery

| Script | Figure |
|---|---|
| `CW-ODMR.py` | ![CW-ODMR spectrum](Figures/CW-ODMR.png) |
| `Ramsey.py` | ![Ramsey fringes](Figures/Ramsey.png) |
| `GSLAC ESLAC con E.py` | ![ESLAC shift with E-field](Figures/ESLAC-GSLAC%20con%20E.png) |
| `AC-Efield.py` | ![AC E-field dynamic response](Figures/Fixed%20bias%20B%20AC%20E.png) |

---

## References

- [Doherty et al., *Physics Reports* 528, 1 (2013)](https://doi.org/10.1016/j.physrep.2013.02.001) — Comprehensive review of NV centre physics.
- [Johansson et al., *Comput. Phys. Commun.* 183, 1760 (2012)](https://doi.org/10.1016/j.cpc.2012.02.021) — QuTiP: Quantum Toolbox in Python.
- Sun, Q. et al. - PySCF: the Python-based simulations of chemistry framework.
- [`NV_E_G_SLAC_explanation.pdf`](7 level Model + Lindbald/NV_E_G_SLAC_explanation.pdf) — Supplementary derivations for the GSLAC/ESLAC regime.

---

## License

Se pone privado en breves...