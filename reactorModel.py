import numpy as np
import sys
import matplotlib.pyplot as plt

import molarMass as mM
import halfLife as hL
import crossSection as cS


def reactorModel(fuelCompo, FPCompo, t_final, n_th_init, n_fa_init, mTot):

    """
    Computes the evolution of the species in the reactor, and compute the power evolution in time
    :param fuelCompo: class
    :param FPCompo: class
            FPCompo.Xe135: percentage of FP that are considered as poison (in this project, we assume only
                Xe135 is produced as a poison)
            FPCompo.FP: percentage of FP that are considered "others fission products" of U235 (assume only one isotope
                different than Xe135 is produced by the fission)
    :param t_final: double
            final time of the simulation in seconds [s]
    :param n_th_init: double
            initial number of thermal neutrons
    :param n_fa_init: double
            initial number of fast neutrons
    :param mTot: double
            total number of kg of fuel at the initial state in [kg]
    :return fuelBurnup: array
            depending on the fuel, returns the burnup
    """

# ===================== Global constants & numerics ========================
    V_core = 10.0                # Reactor core volume [m^3]
    dt = 1e-4                    # Time step [s]
    n_steps = int(t_final / dt) + 1
    time = np.linspace(0.0, t_final, n_steps)

    # Physical constants
    NA = 6.02214076e23           # Avogadro [mol^-1]
    e_charge = 1.602176634e-19   # [C] -> [J/eV]
    m_n = 1.67492749804e-27      # neutron mass [kg]

    # Neutron energies [eV]
    E_th_eV = 0.025              # thermal
    E_fast_eV = 1.0e6            # fast (~1 MeV)

    # Neutron speeds [m/s]  (v = sqrt(2 E / m))
    v_th = np.sqrt(2.0 * E_th_eV * e_charge / m_n)
    v_fast = np.sqrt(2.0 * E_fast_eV * e_charge / m_n)

    # Conversion from microscopic σ [barn] to effective rate coefficient
    # rate ≈ σ [m²] * v [m/s] * n/V [1/m³] * N [-]
    k_th = v_th / V_core * 1e-28
    k_fast = v_fast / V_core * 1e-28  # (non utilisé pour les actinides dans ce modèle simplifié)

    # ======================== Initial fuel inventories ========================
    # Fuel mass fractions in %
    frac_U235 = fuelCompo.U235 / 100.0
    frac_U238 = fuelCompo.U238 / 100.0
    frac_Pu239 = fuelCompo.Pu239 / 100.0
    frac_Th232 = fuelCompo.Th232 / 100.0

    # Numbers of heavy nuclei at t=0
    N_U235 = mTot * frac_U235 / mM.molarMass('U235') * NA
    N_U238 = mTot * frac_U238 / mM.molarMass('U238') * NA
    N_Pu239 = mTot * frac_Pu239 / mM.molarMass('Pu239') * NA
    N_Th232 = mTot * frac_Th232 / mM.molarMass('Th232') * NA  # non utilisé dans ce modèle simple

    # Fission products: Xe-135 (poison) + lumped FP (delayed neutron precursors)
    N_Xe = 0.0
    N_FP = 0.0

    # Neutrons (fast & thermal)
    n_th = n_th_init
    n_fa = n_fa_init

    # =============== Fission product composition & delayed neutrons ===========
    frac_xe = FPCompo.Xe135 / (FPCompo.Xe135 + FPCompo.FP)
    frac_fp = FPCompo.FP / (FPCompo.Xe135 + FPCompo.FP)

    # Prompt neutrons per fission
    nu_prompt = 2.0

    # Effective delayed fraction beta (order-of-magnitude)
    beta_eff = 0.0065
    nu_delayed = beta_eff / (1.0 - beta_eff) * nu_prompt  # delayed neutrons per fission

    # Each fission produces 2 unstable FP, fraction frac_fp are "lumped FP"
    lumps_FP_per_fission = 2.0 * frac_fp
    if lumps_FP_per_fission > 0.0:
        eta_delayed = nu_delayed / lumps_FP_per_fission   # delayed n per decaying FP
    else:
        eta_delayed = 0.0

    # ================= Slowdown & radioactive decay constants =================
    # Fast -> thermal half-time
    lambda_sd = np.log(2.0) / (5e-4)   # [1/s]

    # Lumped FP half-life ~ 1 s
    lambda_FP = np.log(2.0) / 1.0      # [1/s]

    # Xe-135 beta- decay constant
    hl_Xe = hL.halfLife('Xe135', 'BetaMinus')
    if np.isinf(hl_Xe) or hl_Xe <= 0.0:
        lambda_Xe = 0.0
    else:
        lambda_Xe = np.log(2.0) / hl_Xe

    # ==================== Microscopic cross sections (thermal) =================
    E_th = E_th_eV
    sigma_f_U235_th = cS.crossSection('U235', 'Fission', E_th)
    sigma_c_U235_th = cS.crossSection('U235', 'Capture', E_th)
    sigma_f_U238_th = cS.crossSection('U238', 'Fission', E_th)
    sigma_c_U238_th = cS.crossSection('U238', 'Capture', E_th)
    sigma_f_Pu239_th = cS.crossSection('Pu239', 'Fission', E_th)
    sigma_c_Pu239_th = cS.crossSection('Pu239', 'Capture', E_th)
    sigma_c_Xe_th = cS.crossSection('Xe135', 'Capture', E_th)

    # ======================= Leak & control-rod parameters =====================
    # Leak terms (constant)
    SIGMA_LEAK_FAST = 0.5     # [1/s]
    SIGMA_LEAK_TH = 0.05      # [1/s]

    # Control rods: additional absorption, both thermal and fast, in [0;20] 1/s
    SIGMA_CTRL_MAX = 20.0     # [1/s]
    sigma_ctrl_th = 0.0       # thermal control Σ_th
    sigma_ctrl_fast = 0.0     # fast   control Σ_fast

    # ====================== Energies per nuclear event ========================
    # Typical values :  ~200 MeV per fission, split into:
    #  - prompt fission energy (~180 MeV)
    #  - delayed FP decay energy (~10 MeV)
    #  - neutron slowing-down (~1 MeV)
    E_FISSION_PROMPT = 180.0e6 * e_charge  # [J] per fission
    E_FP_DECAY = 10.0e6 * e_charge         # [J] per FP decay
    E_SLOWDOWN = E_fast_eV * e_charge      # [J] per neutron slowed to thermal

    # ============================ Storage arrays ==============================
    burnup = np.zeros(n_steps)       # [MWd/kg]
    P_total_arr = np.zeros(n_steps)  # [W]
    n_th_arr = np.zeros(n_steps)
    n_fa_arr = np.zeros(n_steps)
    N_U235_arr = np.zeros(n_steps)
    N_U238_arr = np.zeros(n_steps)
    N_Pu239_arr = np.zeros(n_steps)
    N_Xe_arr = np.zeros(n_steps)
    sigma_ctrl_th_arr = np.zeros(n_steps)
    sigma_ctrl_fast_arr = np.zeros(n_steps)

    E_cum = 0.0      # cumulative energy [J]

    # ---------- Power setpoint (steady power you want) ----------
    # À ajuster si besoin. Pour toi, un ordre de grandeur raisonnable
    # vu ton pic ~1.5e10 W est P_set ≈ 1e9 W.
    P_set = 1e9      # [W]

    # ---------- P-controller gains for control rods ----------
    # On ne multiplie PLUS par dt dans la loi de contrôle (voir plus bas),
    # donc les gains doivent être petits.
    Kp_th = 1e-11    # gain pour Σ_ctrl_th
    Kp_fast = 5e-12  # gain pour Σ_ctrl_fast

    # ============================ Time integration ============================
    for k in range(n_steps):
        t = k * dt

        # --- store current state ---
        n_th_arr[k] = n_th
        n_fa_arr[k] = n_fa
        N_U235_arr[k] = N_U235
        N_U238_arr[k] = N_U238
        N_Pu239_arr[k] = N_Pu239
        N_Xe_arr[k] = N_Xe
        sigma_ctrl_th_arr[k] = sigma_ctrl_th
        sigma_ctrl_fast_arr[k] = sigma_ctrl_fast

        # ---------------------- Reaction rates [1/s] -------------------------
        # Fissions (thermal)
        Rf_U235 = k_th * sigma_f_U235_th * n_th * N_U235
        Rf_U238 = k_th * sigma_f_U238_th * n_th * N_U238
        Rf_Pu239 = k_th * sigma_f_Pu239_th * n_th * N_Pu239

        # Captures (thermal)
        Rc_U235 = k_th * sigma_c_U235_th * n_th * N_U235
        Rc_U238 = k_th * sigma_c_U238_th * n_th * N_U238
        Rc_Pu239 = k_th * sigma_c_Pu239_th * n_th * N_Pu239

        # Xe-135 capture
        Rcap_Xe = k_th * sigma_c_Xe_th * n_th * N_Xe

        # Total fission rate
        Rf_total = Rf_U235 + Rf_U238 + Rf_Pu239

        # --------------------- Neutron sources & losses ----------------------
        # Fast -> thermal slowdown
        S_slowdown = lambda_sd * n_fa

        # Delayed neutrons from FP decays (assumed directly thermal)
        S_delayed = eta_delayed * lambda_FP * N_FP

        # Total losses of thermal neutrons by nuclear interactions + Xe + leak + rods
        Loss_th_nuc = (Rf_total + Rc_U235 + Rc_U238 + Rc_Pu239 + Rcap_Xe)
        Loss_th_add = (SIGMA_LEAK_TH + sigma_ctrl_th) * n_th

        # ODEs for neutrons
        dn_th_dt = S_slowdown + S_delayed - Loss_th_nuc - Loss_th_add
        dn_fa_dt = (nu_prompt * Rf_total
                    - lambda_sd * n_fa
                    - (SIGMA_LEAK_FAST + sigma_ctrl_fast) * n_fa)

        # ODEs for heavy nuclides (fission + capture)
        dN_U235_dt = -(Rf_U235 + Rc_U235)
        dN_U238_dt = -(Rf_U238 + Rc_U238)
        dN_Pu239_dt = -(Rf_Pu239 + Rc_Pu239)
        dN_Th232_dt = 0.0  # chaîne thorium négligée ici

        # ODEs for fission products
        dN_Xe_dt = 2.0 * frac_xe * Rf_total - lambda_Xe * N_Xe - Rcap_Xe
        dN_FP_dt = 2.0 * frac_fp * Rf_total - lambda_FP * N_FP

        # ----------------------------- Power ---------------------------------
        P_fission = Rf_total * E_FISSION_PROMPT
        P_FP = (lambda_Xe * N_Xe + lambda_FP * N_FP) * E_FP_DECAY
        P_slow = lambda_sd * n_fa * E_SLOWDOWN
        P_total = P_fission + P_FP + P_slow
        P_total_arr[k] = P_total

                # ----------------------- Control rods (feedback) ---------------------
        # On veut P_total ≈ P_set.
        error = P_total - P_set  # >0 : trop de puissance → on insère des barres

        # Correcteur proportionnel simple (PAS de *dt ici)
        sigma_ctrl_th += Kp_th * error
        sigma_ctrl_fast += Kp_fast * error

        # Saturation Σ ∈ [0 ; SIGMA_CTRL_MAX]
        sigma_ctrl_th = max(0.0, min(SIGMA_CTRL_MAX, sigma_ctrl_th))
        sigma_ctrl_fast = max(0.0, min(SIGMA_CTRL_MAX, sigma_ctrl_fast))


        # ---------------------- Euler explicit update ------------------------
        n_th = max(0.0, n_th + dn_th_dt * dt)
        n_fa = max(0.0, n_fa + dn_fa_dt * dt)
        N_U235 = max(0.0, N_U235 + dN_U235_dt * dt)
        N_U238 = max(0.0, N_U238 + dN_U238_dt * dt)
        N_Pu239 = max(0.0, N_Pu239 + dN_Pu239_dt * dt)
        N_Th232 = max(0.0, N_Th232 + dN_Th232_dt * dt)
        N_Xe = max(0.0, N_Xe + dN_Xe_dt * dt)
        N_FP = max(0.0, N_FP + dN_FP_dt * dt)

        # ------------------- Cumulative energy & burnup ----------------------
        E_cum += P_total * dt  # [J]
        burnup[k] = E_cum / (mTot * 1e6 * 86400.0)  # [MWd/kg]

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
    }

    return results


class Fuel:
    def __init__(self):
        # Mass fractions in %
        self.U235 = 3
        self.U238 = 97
        self.Pu239 = 0
        self.Th232 = 0


class FP:
    def __init__(self):
        # Fraction of fission products (in %) that are Xe vs "other FP"
        self.Xe135 = 5
        self.FP = 95


# Example usage + plots
if __name__ == "__main__":
    fuelCompo = Fuel()
    FPCompo = FP()

    # Attention: t_final=200 s avec dt=1e-4 → 2e6 pas → long et lourd.
    # Pour tester, commence par t_final=1 ou 5 s.
    t_final = 200
    results = reactorModel(
        fuelCompo=fuelCompo,
        FPCompo=FPCompo,
        t_final=t_final,
        n_th_init=1e10,
        n_fa_init=0.0,
        mTot=25.0,
    )

    time = results["time"]

    # --------- Figure 1: puissance du réacteur ---------
    plt.figure()
    plt.plot(time, results["P_total"])
    plt.xlabel("Time [s]")
    plt.ylabel("Power [W]")
    plt.title("Reactor power vs time")
    plt.grid(True)

    # --------- Figure 2: neutrons thermiques / rapides ---------
    plt.figure()
    plt.semilogy(time, results["n_th"], label="n_th (thermal)")
    plt.semilogy(time, results["n_fa"], label="n_fast (fast)")
    plt.xlabel("Time [s]")
    plt.ylabel("Neutron population [-]")
    plt.title("Neutron populations vs time")
    plt.legend()
    plt.grid(True)

    # --------- Figure 3: U-235 et Pu-239 ---------
    plt.figure()
    plt.plot(time, results["N_U235"], label="U-235")
    plt.plot(time, results["N_Pu239"], label="Pu-239")
    plt.xlabel("Time [s]")
    plt.ylabel("Number of nuclei [-]")
    plt.title("Fuel evolution (U-235 burnup & Pu-239 buildup)")
    plt.legend()
    plt.grid(True)

    # --------- Figure 4: Xe-135 et contrôles ---------
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.plot(time, results["N_Xe"])
    plt.ylabel("N_Xe [-]")
    plt.title("Xenon-135 build-up")

    plt.subplot(2, 1, 2)
    plt.plot(time, results["sigma_ctrl_th"], label="Σ_ctrl_th")
    plt.plot(time, results["sigma_ctrl_fast"], label="Σ_ctrl_fast")
    plt.xlabel("Time [s]")
    plt.ylabel("Control absorption [1/s]")
    plt.legend()
    plt.grid(True)

    # --------- Figure 5: Burnup ---------
    plt.figure()
    plt.plot(time, results["burnup"])
    plt.xlabel("Time [s]")
    plt.ylabel("Burnup [MWd/kg]")
    plt.title("Fuel burnup vs time")
    plt.grid(True)

    plt.tight_layout()
    plt.show()