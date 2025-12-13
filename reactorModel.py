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
    :return: dict
            time histories (power, neutrons, fuel, Xe, control, burnup, k_eff)
    """
    V_core = 10.0                
    dt = 1e-4                    
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
    k_fast = v_fast / V_core * 1e-28  # ready if you want fast reactions later

    # ======================== Initial fuel inventories ========================
    frac_U235 = fuelCompo.U235 / 100.0
    frac_U238 = fuelCompo.U238 / 100.0
    frac_Pu239 = fuelCompo.Pu239 / 100.0
    frac_Th232 = fuelCompo.Th232 / 100.0

    N_U235 = mTot * frac_U235 / mM.molarMass('U235') * NA
    N_U238 = mTot * frac_U238 / mM.molarMass('U238') * NA
    N_Pu239 = mTot * frac_Pu239 / mM.molarMass('Pu239') * NA
    N_Th232 = mTot * frac_Th232 / mM.molarMass('Th232') * NA  # unused in simplified model

    # Fission products
    N_Xe = 0.0
    N_FP = 0.0

    # Neutrons
    n_th = n_th_init
    n_fa = n_fa_init

    # =============== Fission product composition & delayed neutrons ===========
    frac_xe = FPCompo.Xe135 / (FPCompo.Xe135 + FPCompo.FP)
    frac_fp = FPCompo.FP / (FPCompo.Xe135 + FPCompo.FP)

    nu_prompt = 2.0                # prompt neutrons per fission

    beta_eff = 0.0065              # effective delayed fraction
    nu_delayed = beta_eff / (1.0 - beta_eff) * nu_prompt

    lumps_FP_per_fission = 2.0 * frac_fp
    if lumps_FP_per_fission > 0.0:
        eta_delayed = nu_delayed / lumps_FP_per_fission   # delayed neutrons per decaying FP
    else:
        eta_delayed = 0.0

    # ================= Slowdown & radioactive decay constants =================
    lambda_sd = np.log(2.0) / (5e-4)   # [1/s] fast -> thermal (given)

    lambda_FP = np.log(2.0) / 1.0      # [1/s] FP lumped (T1/2 ~ 1s, hint)

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

    # ======================= Leak & control-rod parameters ====================
    SIGMA_LEAK_FAST = 0.5     # [1/s]
    SIGMA_LEAK_TH = 0.05      # [1/s]

    SIGMA_CTRL_MAX = 20.0    # [1/s]

    # Position de référence des barres : un peu plus absorbant au départ
    # => k_eff(t=0) < 1 (sous-critique, barres plus enfoncées)
    SIGMA_EQ_TH = 14.9       # [1/s]
    SIGMA_EQ_FAST = 14.9     # [1/s]

    # On initialise les Σ_ctrl à leur valeur de référence
    sigma_ctrl_th = SIGMA_EQ_TH
    sigma_ctrl_fast = SIGMA_EQ_FAST

    # ====================== Energies per nuclear event ========================
    E_FISSION_PROMPT = 180.0e6 * e_charge  # [J] per fission
    E_FP_DECAY = 10.0e6 * e_charge         # [J] per FP decay
    E_SLOWDOWN = E_fast_eV * e_charge      # [J] per neutron slowed to thermal

    # ============================ Storage arrays ==============================
    burnup = np.zeros(n_steps)
    P_total_arr = np.zeros(n_steps)
    n_th_arr = np.zeros(n_steps)
    n_fa_arr = np.zeros(n_steps)
    N_U235_arr = np.zeros(n_steps)
    N_U238_arr = np.zeros(n_steps)
    N_Pu239_arr = np.zeros(n_steps)
    N_Xe_arr = np.zeros(n_steps)
    sigma_ctrl_th_arr = np.zeros(n_steps)
    sigma_ctrl_fast_arr = np.zeros(n_steps)
    k_eff_arr = np.zeros(n_steps)          # <-- k_eff(t)

    E_cum = 0.0  # cumulative energy [J]

    # ---------- Power setpoint (with ramp) ----------
    P_nominal = 1e9      # [W] puissance nominale visée
    t_ramp = 50.0        # [s] durée de la rampe 0 -> P_nominal

    # ---------- P-controller (statique autour d'un Σ_eq) ----------
    # Gains assez doux pour éviter l'overshoot
    Kp_th_rel = 3.0      # gain sans dimension pour les thermiques
    Kp_fast_rel = 1.5    # gain sans dimension pour les rapides

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
        Rf_U235 = k_th * sigma_f_U235_th * n_th * N_U235
        Rf_U238 = k_th * sigma_f_U238_th * n_th * N_U238
        Rf_Pu239 = k_th * sigma_f_Pu239_th * n_th * N_Pu239

        Rc_U235 = k_th * sigma_c_U235_th * n_th * N_U235
        Rc_U238 = k_th * sigma_c_U238_th * n_th * N_U238
        Rc_Pu239 = k_th * sigma_c_Pu239_th * n_th * N_Pu239

        Rcap_Xe = k_th * sigma_c_Xe_th * n_th * N_Xe

        Rf_total = Rf_U235 + Rf_U238 + Rf_Pu239

        # --------------------- Neutron sources & losses ----------------------
        S_slowdown = lambda_sd * n_fa
        S_delayed = eta_delayed * lambda_FP * N_FP

        Loss_th_nuc = (Rf_total + Rc_U235 + Rc_U238 + Rc_Pu239 + Rcap_Xe)
        Loss_th_add = (SIGMA_LEAK_TH + sigma_ctrl_th) * n_th

        dn_th_dt = S_slowdown + S_delayed - Loss_th_nuc - Loss_th_add
        dn_fa_dt = (nu_prompt * Rf_total
                    - lambda_sd * n_fa
                    - (SIGMA_LEAK_FAST + sigma_ctrl_fast) * n_fa)

        dN_U235_dt = -(Rf_U235 + Rc_U235)
        dN_U238_dt = -(Rf_U238 + Rc_U238)
        dN_Pu239_dt = -(Rf_Pu239 + Rc_Pu239)
        dN_Th232_dt = 0.0

        dN_Xe_dt = 2.0 * frac_xe * Rf_total - lambda_Xe * N_Xe - Rcap_Xe
        dN_FP_dt = 2.0 * frac_fp * Rf_total - lambda_FP * N_FP

        # ----------------------------- Power ---------------------------------
        P_fission = Rf_total * E_FISSION_PROMPT
        P_FP = (lambda_Xe * N_Xe + lambda_FP * N_FP) * E_FP_DECAY
        P_slow = lambda_sd * n_fa * E_SLOWDOWN
        P_total = P_fission + P_FP + P_slow
        P_total_arr[k] = P_total

        # ----------------------- k_eff instantané ----------------------------
        # Production de neutrons : prompts + retardés
        prod_neutrons = nu_prompt * Rf_total + S_delayed

        # Pertes : absorptions nucléaires + fuites + barres
        loss_leak_ctrl = (SIGMA_LEAK_FAST + sigma_ctrl_fast) * n_fa \
                         + (SIGMA_LEAK_TH + sigma_ctrl_th) * n_th
        loss_abs_nuc = Loss_th_nuc
        loss_total = loss_abs_nuc + loss_leak_ctrl

        if loss_total > 0.0:
            k_eff_inst = prod_neutrons / loss_total
        else:
            k_eff_inst = 0.0

        k_eff_arr[k] = k_eff_inst

        # ----------------------- Control rods (feedback) ---------------------
        # Consigne de puissance rampée 0 -> P_nominal sur t_ramp
        if t < t_ramp:
            P_set_t = P_nominal * (t / t_ramp)
        else:
            P_set_t = P_nominal

        # Erreur relative de puissance (normalisée par P_nominal)
        if P_nominal > 0.0:
            rel_error = (P_total - P_set_t) / P_nominal
        else:
            rel_error = 0.0

        # Loi de commande "statique" autour de Σ_eq :
        sigma_ctrl_th = SIGMA_EQ_TH + Kp_th_rel * rel_error
        sigma_ctrl_fast = SIGMA_EQ_FAST + Kp_fast_rel * rel_error

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
        "k_eff": k_eff_arr,          # <-- nouveau
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

    t_final = 100.0
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

    # --------- Figure 2: k_eff vs time ---------
    plt.figure()
    plt.plot(time, results["k_eff"])
    plt.xlabel("Time [s]")
    plt.ylabel("k_eff [-]")
    plt.title("Effective multiplication factor vs time")
    plt.grid(True)
    # si tu veux zoomer autour de 1 :
    # plt.ylim(0.9, 1.1)

    plt.show()
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