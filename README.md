
# LMECA2600 – Project 2025  
## Simplified thermal PWR model (Python code)

This project implements a **time–dependent model of a simplified Pressurized Water Reactor (PWR)** with two neutron energy groups (fast / thermal), basic fuel depletion and xenon poisoning.  

The goal of the main function `reactorModel` is to simulate the short–time behaviour of the reactor and to **reach and maintain a nominal power of ≈ 1 GW** by acting on the control rods, modelled as variable macroscopic absorption cross sections.

Time integration is done with an explicit Euler scheme and a constant time step  
\(\Delta t = 10^{-4}\,\text{s}\), as required by the project statement.

---

## File overview

### 1. `crossSection.py`

**Main function**

```python
crossSection(X, Transfo, E_neutron)
```

**Goal**

Returns the **microscopic cross section** σ(X, Transfo, E) in **barn** for a given nuclide and reaction type, at the requested neutron energy(ies).

**Arguments**

- `X` (str): nuclide symbol, e.g. `"U235"`, `"U238"`, `"Pu239"`, `"Xe135"`, `"n"`, …  
- `Transfo` (str): reaction type, in  
  `["Fission", "Capture"]`.  
- `E_neutron` (float or array-like): neutron energy in **eV**, within the project range  
  `[1e-5; 2e7]` eV.

**Behaviour / implementation**

- Internally, the file contains a **small database** of tabulated cross sections (taken from ENDF/B data, simplified) for the nuclides needed in the model, separated for **fission** and **radiative capture**.
- For each nuclide and reaction, the code performs a simple interpolation (or returns a constant representative value) between a **thermal–like energy** (~0.025 eV) and a **fast–like energy** (~1 MeV), according to the needs of the two–group model.
- If a requested nuclide or transformation is not implemented, the function prints a warning and returns 0.0.

**Additional helpers**

- Internal helper functions (e.g. `_get_table(X, Transfo)` or `_interpolate(E, E_grid, sigma_grid)`) may be defined for clarity.  
  They are only used inside `crossSection` and are not called directly in the tests.

---

### 2. `halfLife.py`

**Main function**

```python
halfLife(X, Transfo)
```

**Goal**

Returns the **physical half-life** \(T_{1/2}\) of a nuclide for a specific decay mode, in **seconds**.

**Arguments**

- `X` (str): nuclide symbol, e.g. `"U239"`, `"Np239"`, `"Pu239"`, `"Xe135"`, …  
- `Transfo` (str): decay mode in  
  `["Alpha", "BetaMinus", "BetaPlus", "Gamma"]`.

**Behaviour / implementation**

- The file contains a small **hard-coded database** extracted from the JAEA Chart of the Nuclides (2014) for all nuclides that appear in the kinetic scheme (U, Np, Pu, Xe).  
  Times are converted to **seconds** (H, M, y → s).
- For short–lived isomeric states that only serve as intermediate steps, very small half-lives (or `np.inf` for effectively stable states) are used, consistent with the project simplifications.
- If `X` or `Transfo` is not present in the database, the function prints a warning and returns `np.inf`.

---

### 3. `molarMass.py`

**Main function**

```python
molarMass(X)
```

**Goal**

Returns the **molar mass** \(M_X\) of nuclide or nucleon `X` in **kg/mol**.

**Arguments**

- `X` (str): nuclide symbol or `"n"` for a free neutron, for example `"U235"`, `"U238"`, `"Pu239"`, `"Th232"`, `"Xe135"`, `"n"`.

**Behaviour / implementation**

- The molar masses (in g/mol) are taken from the JAEA Nuclide Chart database and converted to kg/mol.
- Only the isotopes actually used in the model are stored (Th-232, U-233/235/238, U-239, Np-239, Pu-239/240, Xe-135, neutron).
- If `X` is not implemented, the function prints a warning and returns 0.0.

---

### 4. `reactorModel.py`

**Main function**

```python
reactorModel(fuelCompo, FPCompo, t_final, n_th_init, n_fa_init, mTot)
```

This is the **core of the project** and the file that will be used for the assessment.

**Goal**

Simulate the **time evolution** of a simplified PWR over short times (a few hundred seconds), including:

- two neutron energy groups (fast / thermal),
- fission of the main fuel isotopes,
- production and decay of lumped fission products (FP) and xenon-135,
- conversion of fast neutrons to thermal neutrons,
- fuel depletion and xenon poisoning,
- neutron leaks and **control rods** modelled as variable macroscopic capture cross sections  
  \(\Sigma_\text{th}, \Sigma_\text{fast} \in [0; 20]\ \text{s}^{-1}\),
- calculation of reactor **power** and **burnup**,
- simple feedback law on the control rods in order to reach and maintain a **nominal power of ~1 GW**.

**Arguments**

- `fuelCompo` (class `Fuel`): initial fuel composition in mass percentages (U-235, U-238, Pu-239, Th-232).  
- `FPCompo` (class `FP`): partition of fission products between Xe-135 and “other FP”:
  - `FPCompo.Xe135` – percentage of FP considered as xenon,
  - `FPCompo.FP` – percentage of FP considered as “other isotopes” (single lumped nuclide, delayed neutron precursors).
- `t_final` (float): total simulation time [s].  
- `n_th_init` (float): initial number of thermal neutrons [-].  
- `n_fa_init` (float): initial number of fast neutrons [-].  
- `mTot` (float): total initial fuel mass [kg] (25 kg as required in the statement).

**Main modelling choices**

- Two-group model with **E_fast ≈ 1 MeV** and **E_th ≈ 0.025 eV**, and a half-time of 5×10⁻⁴ s for fast neutrons to thermalize.
- Each fission produces:
  - ν = 2 prompt neutrons,  
  - two unstable FP with a half-life of ≈ 1 s (one branch Xe-135, the other lumped “FP”).  
- “Other FP” are the **only source of delayed neutrons** (lumped precursor group).
- Three energy contributions are included in the power:
  1. fission energy (≈ 180 MeV per fission),
  2. decay/stabilization of fission products (≈ 10 MeV per FP),
  3. energy released when fast neutrons are slowed down (≈ E_fast per neutron).
- Control rods and leakage are represented by **macroscopic absorption cross sections** added to the loss term of each neutron group:

  \[
  \frac{dn_\text{th}}{dt} = \dots - \big(\Sigma_\text{leak,th} + \Sigma_\text{ctrl,th}\big)\,n_\text{th},\qquad
  \frac{dn_\text{fast}}{dt} = \dots - \big(\Sigma_\text{leak,fast} + \Sigma_\text{ctrl,fast}\big)\,n_\text{fast}
  \]

  with \(\Sigma_\text{ctrl,th},\Sigma_\text{ctrl,fast} \in [0; 20]\ \mathrm{s}^{-1}\) as required.

- A **simple proportional feedback** on power adjusts the control macroscopic cross sections:

  - a power setpoint `P_set(t)` is ramped from 0 to 1 GW over a given time `t_ramp`,
  - at each time step, the relative error `(P - P_set)/P_nominal` is used to slightly increase or decrease `Σ_ctrl_th` and `Σ_ctrl_fast`, within [0; 20] s⁻¹,
  - this brings the reactor towards a **quasi-steady power close to 1 GW** (small overshoot allowed).

- The effective multiplication factor \(k_\text{eff}(t)\) is estimated as

  \[
  k_\text{eff} = \frac{\text{neutrons produced by fission and delayed sources}}
                      {\text{all absorptions + leaks + control rods}}.
  \]

**Outputs**

`reactorModel` returns a dictionary:

```python
results = {
    "time": time,
    "burnup": burnup,
    "P_total": P_total_arr,
    "n_th": n_th_arr,
    "n_fa": n_fa_arr,
    "N_U235": N_U235_arr,
    "N_U238": N_U238_arr,
    "N_Pu239": N_Pu239_arr,
    "N_Xe": N_Xe_arr,
    "sigma_ctrl_th": sigma_ctrl_th_arr,
    "sigma_ctrl_fast": sigma_ctrl_fast_arr,
    "k_eff": k_eff_arr,
}
```

**Helper classes**

- `class Fuel`: container for initial fuel mass fractions (**U-235, U-238, Pu-239, Th-232**).  
- `class FP`: container for FP partition between **Xe135** and **other FP**.

---

### 5. `project.py` (optional driver script)

This script is used to perform **several analyses** based on the `reactorModel` function, for example:

- **Base case**: nominal conditions (25 kg fuel, standard 3% enrichment, default control parameters) and plots of power, k_eff, Xe-135, burnup.
- **Sensitivity to control parameters**: change the proportional gain, the ramp time, or the initial control cross sections and compare the transient and the overshoot.
- **Effect of xenon poisoning**: compare runs with and without Xe-135 absorption (σ_Xe = 0) to show the impact on power and control rod position.
- **Effect of fuel composition**: vary Pu-239 content or enrichment and see how much control-rod worth is required to keep 1 GW.

A typical structure is:

```python
from reactorModel import reactorModel, Fuel, FP
import matplotlib.pyplot as plt

def run_base_case():
    fuel = Fuel()
    fp = FP()
    res = reactorModel(fuel, fp,
                       t_final=100.0,
                       n_th_init=1e10,
                       n_fa_init=0.0,
                       mTot=25.0)
    # plotting commands …

def main():
    run_base_case()
    # possibly call other scenarios

if __name__ == "__main__":
    main()
```

---

## How to run

- To **test individual modules**, you can run each file directly:

```bash
python crossSection.py
python halfLife.py
python molarMass.py
python reactorModel.py
```

Each file contains a small `if __name__ == "__main__":` block that runs a simple example.

- To reproduce all the analyses used in your presentation, run:

```bash
python project.py
```

which will call `reactorModel` with the different parameter sets and generate the corresponding plots.

---

## Summary of main assumptions

- Two-group neutron model (fast / thermal) with fixed representative energies.
- Only a limited set of nuclides are explicitly tracked: U-235, U-238, Pu-239, Th-232, Xe-135 and one lumped “other FP” precursor.
- 2 FP per fission, one branch Xe-135, one branch “other FP” (delayed neutron precursor).
- Single group of delayed neutrons, effective fraction β_eff taken from literature.
- Control rods and leaks represented by macroscopic capture cross sections Σ_th, Σ_fast in [0; 20] s⁻¹.
- Explicit Euler time integration with Δt = 10⁻⁴ s.
- Cross sections and molar masses are simplified, representative values taken from nuclear data libraries.
