import numpy as np
import matplotlib.pyplot as plt

from reactorModel import reactorModel, Fuel, FP


def find_optimal_sigma_thermal(P_target,
                               t_scan=5.0,
                               mTot=25.0,
                               n_th_init=1e10,
                               n_fa_init=1e10):
    """
    Balaye plusieurs valeurs de Sigma_thermal pour approcher une puissance cible P_target.
    Retourne (Sigma_thermal_opt, resultats_simulation).
    """

    fuelCompo = Fuel()
    FPCompo = FP()

    # On balaie Sigma_thermal entre 0 et 20 [1/s]
    sigma_candidates = np.linspace(0.0, 20.0, 21)

    best_sigma = None
    best_diff = None
    best_res = None

    for sigma_th in sigma_candidates:
        res = reactorModel(
            fuelCompo=fuelCompo,
            FPCompo=FPCompo,
            t_final=t_scan,
            n_th_init=n_th_init,
            n_fa_init=n_fa_init,
            mTot=mTot,
            Sigma_fast=0.0,
            Sigma_thermal=sigma_th,
        )

        P = res["P"]
        t = res["t"]

        # Moyenne de la puissance sur le dernier cinquième du temps
        n = len(t)
        idx_start = int(0.8 * n)
        P_mean = P[idx_start:].mean()

        diff = abs(P_mean - P_target)

        print(f"Sigma_thermal = {sigma_th:.2f}  ->  P_mean ~ {P_mean:.3e} W  (diff = {diff:.3e})")

        if (best_diff is None) or (diff < best_diff):
            best_diff = diff
            best_sigma = sigma_th
            best_res = res

    print("\n---------------------------------------")
    print(f"Meilleur Sigma_thermal ≈ {best_sigma:.3f} [1/s]")
    print(f"Puissance moyenne correspondante ≈ {best_diff + P_target:.3e} W")
    print("---------------------------------------\n")

    return best_sigma, best_res


def run_final_simulation(P_target=1e8,
                         t_final=20.0,
                         mTot=25.0,
                         n_th_init=1e10,
                         n_fa_init=1e10):
    """
    1) Cherche un Sigma_thermal qui donne une puissance proche de P_target sur 5 s.
    2) Relance une simulation plus longue avec ce Sigma_thermal.
    3) Trace les courbes utiles pour l'analyse.
    """

    # Étape 1 : recherche de Sigma_thermal optimal sur 5 s
    sigma_opt, _ = find_optimal_sigma_thermal(
        P_target=P_target,
        t_scan=5.0,
        mTot=mTot,
        n_th_init=n_th_init,
        n_fa_init=n_fa_init,
    )

    fuelCompo = Fuel()
    FPCompo = FP()

    # Étape 2 : simulation finale plus longue
    res = reactorModel(
        fuelCompo=fuelCompo,
        FPCompo=FPCompo,
        t_final=t_final,
        n_th_init=n_th_init,
        n_fa_init=n_fa_init,
        mTot=mTot,
        Sigma_fast=0.0,
        Sigma_thermal=sigma_opt,
    )

    t = res["t"]
    P = res["P"]
    n_fast = res["n_fast"]
    n_thermal = res["n_thermal"]
    Xe = res["Xe135"]
    U235 = res["U235"]
    Pu239 = res["Pu239"]

    # Normalisations pour les plots (optionnel)
    U2350 = U235[0] if U235[0] > 0 else 1.0
    Pu2390 = Pu239[0] if Pu239[0] > 0 else 1.0

    # --- Figure 1 : puissance ---
    plt.figure()
    plt.plot(t, P)
    plt.xlabel("Temps (s)")
    plt.ylabel("Puissance (W)")
    plt.title(f"Puissance du réacteur (Sigma_th ≈ {sigma_opt:.2f} 1/s)")
    plt.grid(True)

    # --- Figure 2 : neutrons rapides / thermiques ---
    plt.figure()
    plt.plot(t, n_fast, label="n_fast")
    plt.plot(t, n_thermal, label="n_thermal")
    plt.xlabel("Temps (s)")
    plt.ylabel("Nombre de neutrons")
    plt.title("Évolution des neutrons rapides et thermiques")
    plt.legend()
    plt.grid(True)

    # --- Figure 3 : Xe-135 ---
    plt.figure()
    plt.plot(t, Xe)
    plt.xlabel("Temps (s)")
    plt.ylabel("N_Xe135 (arbitraire)")
    plt.title("Évolution du poison Xe-135")
    plt.grid(True)

    # --- Figure 4 : fractions U235 / Pu239 ---
    plt.figure()
    plt.plot(t, U235 / U2350, label="U235 / U235(0)")
    plt.plot(t, Pu239 / Pu2390, label="Pu239 / Pu239(0)")
    plt.xlabel("Temps (s)")
    plt.ylabel("Fraction normalisée")
    plt.title("Évolution de U235 et Pu239 (normalisées)")
    plt.legend()
    plt.grid(True)

    plt.show()


if __name__ == "__main__":
    # Puissance cible réaliste pour ce cœur très sous-critique (~1 W)
    P_target = 1.0   # en Watt
    run_final_simulation(P_target=P_target)

