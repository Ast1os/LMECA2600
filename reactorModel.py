import numpy as np
from crossSection import crossSection
from halfLife import halfLife
from molarMass import molarMass


def reactorModel(fuelCompo, FPCompo, t_final, n_th_init, n_fa_init, mTot,
                 Sigma_fast=0.0, Sigma_thermal=0.0):
    """
    Simulateur simplifié d'un réacteur PWR à deux groupes d'énergie (rapide / thermique).

    Paramètres
    ----------
    fuelCompo : objet
        Doit contenir les attributs (en pourcentage massique) :
        - fuelCompo.U235
        - fuelCompo.U238

    FPCompo : objet
        Doit contenir les attributs :
        - FPCompo.Xe135  (valeur initiale de Xe-135, en nombre d'atomes ou arbitraire)
        - FPCompo.FP     (non utilisé explicitement dans ce modèle simplifié)

    t_final : float
        Temps final de la simulation (s).

    n_th_init : float
        Nombre initial de neutrons thermiques.

    n_fa_init : float
        Nombre initial de neutrons rapides.

    mTot : float
        Masse totale de combustible (kg). L’énoncé recommande 25 kg.

    Sigma_fast : float, optionnel
        Section efficace macroscopique de capture pour les neutrons rapides [1/s],
        représentant fuites + barres de contrôle agissant sur les rapides.

    Sigma_thermal : float, optionnel
        Section efficace macroscopique de capture pour les neutrons thermiques [1/s],
        représentant fuites + barres de contrôle agissant sur les thermiques.

    Retour
    ------
    dict
        Contient l'évolution temporelle des grandeurs suivantes :
        - "t"        : temps
        - "n_fast"   : neutrons rapides
        - "n_thermal": neutrons thermiques
        - "U235"     : population de noyaux U-235
        - "U238"     : population de noyaux U-238
        - "U239"     : population de noyaux U-239
        - "Np239"    : population de noyaux Np-239
        - "Pu239"    : population de noyaux Pu-239
        - "precursors": précurseurs de neutrons retardés (FP autres que Xe)
        - "Xe135"    : population de Xe-135
        - "P"        : puissance totale (W)
    """

    # ------------------ Paramètres physiques globaux ------------------

    V_core = 10.0               # m^3 (volume du cœur)
    dt = 1e-4                   # s (pas de temps, imposé par l'énoncé)
    E_fast = 1e6                # eV  (ordre de grandeur ~1 MeV)
    E_th = 0.025                # eV  (neutrons thermiques)
    nu = 2.0                    # neutrons prompts par fission
    y_Xe = 0.25                 # fraction des FP qui sont du Xe-135 (choix raisonnable)

    # Neutrons retardés : T1/2 ~ 1 s => lambda_prec :
    lambda_prec = np.log(2.0) / 1.0  # s^-1

    # Thermalisation des neutrons rapides : T1/2 ~ 5e-4 s :
    T12_thermalisation = 5e-4  # s
    lambda_th = np.log(2.0) / T12_thermalisation

    # Demi-vie du Xe-135 (utilise halfLife.py, en s) :
    T12_Xe = halfLife("Xe135", "BetaMinus")
    lambda_Xe = np.log(2.0) / T12_Xe

    # Chaîne U238 -> U239 -> Np239 -> Pu239
    T12_U239 = halfLife("U239", "BetaMinus")
    T12_Np239 = halfLife("Np239", "BetaMinus")
    lambda_U239 = np.log(2.0) / T12_U239
    lambda_Np239 = np.log(2.0) / T12_Np239

    # Vitesse moyenne des neutrons (ordre de grandeur)
    v_fast = 1.4e7   # m/s (neutrons MeV)
    v_th = 2.2e3     # m/s (neutrons thermiques)

    # Énergie par fission (on regroupe fission + stabilisation FF + ralentissement)
    # ~ 205 MeV par fission :
    E_per_fission = 205e6 * 1.602e-19  # J

    # ------------------ Pré-calcul des sections efficaces (barn -> m^2) ------------------

    barn_to_m2 = 1e-28

    # U-235
    sig_f_U235_th = crossSection("U235", "Fission", E_th) * barn_to_m2
    sig_f_U235_fa = crossSection("U235", "Fission", E_fast) * barn_to_m2
    sig_c_U235_th = crossSection("U235", "Capture", E_th) * barn_to_m2
    sig_c_U235_fa = crossSection("U235", "Capture", E_fast) * barn_to_m2

    # U-238 (fission sur rapides négligée ici pour simplifier)
    sig_f_U238_th = 0.0
    sig_f_U238_fa = 0.0
    sig_c_U238_th = crossSection("U238", "Capture", E_th) * barn_to_m2
    sig_c_U238_fa = crossSection("U238", "Capture", E_fast) * barn_to_m2

    # Pu-239
    sig_f_Pu239_th = crossSection("Pu239", "Fission", E_th) * barn_to_m2
    sig_f_Pu239_fa = crossSection("Pu239", "Fission", E_fast) * barn_to_m2
    sig_c_Pu239_th = crossSection("Pu239", "Capture", E_th) * barn_to_m2
    sig_c_Pu239_fa = crossSection("Pu239", "Capture", E_fast) * barn_to_m2

    # Xe-135 (poison neutronique)
    sig_c_Xe_th = crossSection("Xe135", "Capture", E_th) * barn_to_m2
    sig_c_Xe_fa = crossSection("Xe135", "Capture", E_fast) * barn_to_m2

    # ------------------ Discrétisation temporelle ------------------

    Nsteps = int(t_final / dt) + 1

    t = np.zeros(Nsteps)
    n_fast = np.zeros(Nsteps)
    n_thermal = np.zeros(Nsteps)
    P = np.zeros(Nsteps)

    # Noyaux
    N_U235 = np.zeros(Nsteps)
    N_U238 = np.zeros(Nsteps)
    N_U239 = np.zeros(Nsteps)
    N_Np239 = np.zeros(Nsteps)
    N_Pu239 = np.zeros(Nsteps)

    N_prec = np.zeros(Nsteps)
    N_Xe = np.zeros(Nsteps)

    # ------------------ Conditions initiales ------------------

    t[0] = 0.0
    n_fast[0] = n_fa_init
    n_thermal[0] = n_th_init

    # Masse molaire (kg/mol)
    MU235 = molarMass("U235")
    MU238 = molarMass("U238")

    # Nombre d'Avogadro
    NA = 6.02214076e23

    # On suppose que le combustible est uniquement U235 + U238
    m_U235 = (fuelCompo.U235 / 100.0) * mTot
    m_U238 = (fuelCompo.U238 / 100.0) * mTot

    N_U235[0] = m_U235 / MU235 * NA
    N_U238[0] = m_U238 / MU238 * NA

    # Pas d'actinides issus du breeding au départ
    N_U239[0] = 0.0
    N_Np239[0] = 0.0
    N_Pu239[0] = 0.0

    # Fission products initiaux
    # Pour simplifier, on peut prendre la valeur de FPCompo.Xe135 comme nombre initial
    N_Xe[0] = FPCompo.Xe135
    N_prec[0] = 0.0

    # ------------------ Boucle temporelle (Euler explicite) ------------------

    for k in range(Nsteps - 1):
        t[k+1] = t[k] + dt

        # ---- Densité de neutrons et flux (0D) ----
        n_fa_k = n_fast[k]
        n_th_k = n_thermal[k]

        # Densité (neutrons / m^3)
        n_fa_dens = n_fa_k / V_core
        n_th_dens = n_th_k / V_core

        # Flux neutronique (n / (m^2 s)) ~ n_dens * v
        phi_fast = n_fa_dens * v_fast
        phi_th = n_th_dens * v_th

        # ---- Populations de noyaux à l'instant k ----
        U235_k = N_U235[k]
        U238_k = N_U238[k]
        U239_k = N_U239[k]
        Np239_k = N_Np239[k]
        Pu239_k = N_Pu239[k]
        Xe_k = N_Xe[k]
        prec_k = N_prec[k]

        # ===================== RÉACTIONS NEUTRONIQUES =====================

        # --- Fissions U-235 ---
        R_f_U235_th = sig_f_U235_th * phi_th * U235_k
        R_f_U235_fa = sig_f_U235_fa * phi_fast * U235_k

        # --- Fissions Pu-239 ---
        R_f_Pu239_th = sig_f_Pu239_th * phi_th * Pu239_k
        R_f_Pu239_fa = sig_f_Pu239_fa * phi_fast * Pu239_k

        # --- (Optionnel) fission U238 négligée dans ce modèle ---
        R_f_U238_th = sig_f_U238_th * phi_th * U238_k
        R_f_U238_fa = sig_f_U238_fa * phi_fast * U238_k

        # Total fissions (pour neutrons & FP)
        R_f_total = (R_f_U235_th + R_f_U235_fa +
                     R_f_Pu239_th + R_f_Pu239_fa +
                     R_f_U238_th + R_f_U238_fa)

        # --- Captures U-235 ---
        R_c_U235_th = sig_c_U235_th * phi_th * U235_k
        R_c_U235_fa = sig_c_U235_fa * phi_fast * U235_k

        # --- Captures U-238 ---
        R_c_U238_th = sig_c_U238_th * phi_th * U238_k
        R_c_U238_fa = sig_c_U238_fa * phi_fast * U238_k

        # --- Captures Pu-239 ---
        R_c_Pu239_th = sig_c_Pu239_th * phi_th * Pu239_k
        R_c_Pu239_fa = sig_c_Pu239_fa * phi_fast * Pu239_k

        # --- Captures Xe-135 ---
        R_c_Xe_th = sig_c_Xe_th * phi_th * Xe_k
        R_c_Xe_fa = sig_c_Xe_fa * phi_fast * Xe_k

        # Pertes par absorption (chaque réaction consomme un neutron)
        Loss_fast_abs = (R_f_U235_fa + R_c_U235_fa +
                         R_f_Pu239_fa + R_c_Pu239_fa +
                         R_f_U238_fa + R_c_U238_fa +
                         R_c_Xe_fa)

        Loss_th_abs = (R_f_U235_th + R_c_U235_th +
                       R_f_Pu239_th + R_c_Pu239_th +
                       R_f_U238_th + R_c_U238_th +
                       R_c_Xe_th)

        # ===================== CINÉTIQUE DES NEUTRONS =====================

        # Neutrons prompts (toutes fissions → neutrons rapides)
        prod_prompts = nu * R_f_total

        # Neutrons retardés (issus des précurseurs)
        prod_delayed = lambda_prec * prec_k

        # Équation pour n_fast
        dn_fast_dt = (
            prod_prompts          # neutrons prompts
            + prod_delayed        # neutrons retardés
            - lambda_th * n_fa_k  # thermalisation vers n_thermal
            - Sigma_fast * n_fa_k # fuites + barres (rapides)
            - Loss_fast_abs       # absorptions sur U/Pu/U238/Xe (rapides)
        )

        # Équation pour n_thermal
        dn_thermal_dt = (
            lambda_th * n_fa_k    # neutrons qui se thermalisent
            - Sigma_thermal * n_th_k  # fuites + barres (thermiques)
            - Loss_th_abs         # absorptions sur U/Pu/U238/Xe (thermiques)
        )

        # ===================== CINÉTIQUE DES NOYAUX =====================

        # U-235 : perte par fission + capture
        dN_U235_dt = -(
            R_f_U235_th + R_f_U235_fa +
            R_c_U235_th + R_c_U235_fa
        )

        # U-238 : perte par capture (fission U-238 négligée ici)
        R_c_U238_total = R_c_U238_th + R_c_U238_fa
        dN_U238_dt = -R_c_U238_total

        # U-239 : produit par capture sur U-238, décroît par beta-
        dN_U239_dt = R_c_U238_total - lambda_U239 * U239_k

        # Np-239 : produit par décroissance de U-239, décroît par beta-
        dN_Np239_dt = lambda_U239 * U239_k - lambda_Np239 * Np239_k

        # Pu-239 : produit par décroissance de Np-239, perte par fission + capture
        dN_Pu239_dt = (
            lambda_Np239 * Np239_k
            - (R_f_Pu239_th + R_f_Pu239_fa + R_c_Pu239_th + R_c_Pu239_fa)
        )

        # Précurseurs (FP autres que Xe-135)
        dN_prec_dt = (
            2.0 * R_f_total * (1.0 - y_Xe)  # 2 FP par fission, fraction non-Xe
            - lambda_prec * prec_k          # décroissance → neutrons retardés
        )

        # Xe-135 : produit par fraction des FP, perte par capture + décroissance
        dN_Xe_dt = (
            2.0 * R_f_total * y_Xe    # 2 FP par fission, fraction Xe
            - (R_c_Xe_th + R_c_Xe_fa) # capture neutronique
            - lambda_Xe * Xe_k        # décroissance beta-
        )

        # ===================== PUISSANCE =====================

        P[k] = R_f_total * E_per_fission  # W

        # ===================== INTÉGRATION EULER =====================

        n_fast[k+1] = max(n_fa_k + dt * dn_fast_dt, 0.0)
        n_thermal[k+1] = max(n_th_k + dt * dn_thermal_dt, 0.0)

        N_U235[k+1] = max(U235_k + dt * dN_U235_dt, 0.0)
        N_U238[k+1] = max(U238_k + dt * dN_U238_dt, 0.0)
        N_U239[k+1] = max(U239_k + dt * dN_U239_dt, 0.0)
        N_Np239[k+1] = max(Np239_k + dt * dN_Np239_dt, 0.0)
        N_Pu239[k+1] = max(Pu239_k + dt * dN_Pu239_dt, 0.0)

        N_prec[k+1] = max(prec_k + dt * dN_prec_dt, 0.0)
        N_Xe[k+1] = max(Xe_k + dt * dN_Xe_dt, 0.0)

    # Dernier point de puissance
    P[-1] = P[-2]

    return {
        "t": t,
        "n_fast": n_fast,
        "n_thermal": n_thermal,
        "U235": N_U235,
        "U238": N_U238,
        "U239": N_U239,
        "Np239": N_Np239,
        "Pu239": N_Pu239,
        "precursors": N_prec,
        "Xe135": N_Xe,
        "P": P,
    }


# Petites classes utilitaires pour tester vite fait le modèle
class Fuel:
    def __init__(self):
        # pourcentages massiques
        self.U235 = 3.0
        self.U238 = 97.0


class FP:
    def __init__(self):
        # valeurs arbitraires pour l'état initial
        self.Xe135 = 0.0
        self.FP = 0.0


if __name__ == "__main__":
    # Exemple de test rapide sans contrôle (Sigma = 0)
    fuelCompo = Fuel()
    FPCompo = FP()
    res = reactorModel(
        fuelCompo=fuelCompo,
        FPCompo=FPCompo,
        t_final=0.01,
        n_th_init=1e10,
        n_fa_init=1e10,
        mTot=25.0,
        Sigma_fast=0.0,
        Sigma_thermal=0.0,
    )
    print("Simulation terminée. Puissance finale (W) :", res["P"][-1])
