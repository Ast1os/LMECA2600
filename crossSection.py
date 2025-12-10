# crossSection.py
# Calcul simplifié des sections efficaces microscopiques en fonction de l'énergie du neutron.
# Les valeurs numériques sont approximatives mais plausibles et suffisantes
# pour un modèle pédagogique à deux groupes d'énergie (rapide / thermique).

import numpy as np


def crossSection(X, Transfo, E_neutron):
    """
    Computes the cross section for a specific energy level
    ENDF database: https://www-nds.iaea.org/exfor/endf.htm

    ----------------
    Parameters
    ----------
    X : str
        Nuclide or nucleon (for a neutron), following the atomic notation of the element.
        Example: for Uranium 235, X = 'U235'.
        Supported nuclides in this simplified model:
        ['Th232', 'Th233', 'Pa233',
         'U233', 'U235', 'U236', 'U237', 'U238', 'U239',
         'Np239',
         'Pu239', 'Pu240',
         'Xe135']

    Transfo : str
        Name of the considered transformation.
        Must be one of: ['Fission', 'Capture'].

    E_neutron : float, list, or numpy.ndarray
        Incident neutron energy (or array of energies) in [eV].
        Valid range: [1e-5 ; 2e7] eV.

    Returns
    -------
    float, list, or numpy.ndarray
        Value(s) of the microscopic cross section in [barn]
        corresponding to the input parameters.
        Le type de sortie suit celui d'entrée :
        - si E_neutron est un scalaire  -> float
        - si E_neutron est une list     -> list de floats
        - si E_neutron est un np.ndarray -> np.ndarray de floats

    Notes
    -----
    Modèle très simplifié :
    - Deux groupes d'énergie :
        * "thermique"  : E < 0.5 eV  (~0.025 eV dans le projet)
        * "rapide"     : E >= 0.5 eV (~1 MeV dans le projet)
    - Dans chaque groupe, la section efficace est supposée constante.
    """

    # ========================== Vérification des arguments ==========================

    allowed_nuclides = [
        'Th232', 'Th233', 'Pa233',
        'U233', 'U235', 'U236', 'U237', 'U238', 'U239',
        'Np239',
        'Pu239', 'Pu240',
        'Xe135'
    ]

    allowed_transfo = ['Fission', 'Capture']

    if X not in allowed_nuclides:
        print(f"\nWARNING : There is no database for element {X}. "
              f"\nPlease check function information\n")

    if Transfo not in allowed_transfo:
        raise ValueError(f"Transformation '{Transfo}' not supported. "
                         f"Use one of {allowed_transfo}.")

    # Garde le type d'entrée pour le retour
    is_scalar = np.isscalar(E_neutron)
    is_list = isinstance(E_neutron, list)

    # Convertit en array numpy pour faciliter les calculs
    E = np.atleast_1d(np.array(E_neutron, dtype=float))

    # Vérifie la plage d'énergie
    E_min_allowed = 1e-5
    E_max_allowed = 2e7
    if (E < E_min_allowed).any() or (E > E_max_allowed).any():
        print("\nWARNING : Some neutron energies are outside the recommended range "
              f"[{E_min_allowed:.1e} ; {E_max_allowed:.1e}] eV.\n")

    # ========================== Base de données simplifiée ==========================

    # Valeurs approximatives (en barn) pour deux groupes d'énergie :
    #   - "th" : thermique  (E < 0.5 eV)
    #   - "fa" : rapide     (E >= 0.5 eV)
    #
    # Pour les réactions physiquement quasi-impossibles, on met 0.0.

    sigma_th = {}  # (X, Transfo) -> sigma_thermal
    sigma_fa = {}  # (X, Transfo) -> sigma_fast

    # --- Thorium / Protactinium / U233 (peu utilisés mais définis pour cohérence) ---
    sigma_th[('Th232', 'Capture')] = 7.0
    sigma_fa[('Th232', 'Capture')] = 0.3
    sigma_th[('Th232', 'Fission')] = 0.0
    sigma_fa[('Th232', 'Fission')] = 0.0

    sigma_th[('Th233', 'Capture')] = 20.0
    sigma_fa[('Th233', 'Capture')] = 0.5
    sigma_th[('Th233', 'Fission')] = 0.0
    sigma_fa[('Th233', 'Fission')] = 0.5

    sigma_th[('Pa233', 'Capture')] = 200.0
    sigma_fa[('Pa233', 'Capture')] = 1.0
    sigma_th[('Pa233', 'Fission')] = 0.0
    sigma_fa[('Pa233', 'Fission')] = 1.5

    sigma_th[('U233', 'Capture')] = 60.0
    sigma_fa[('U233', 'Capture')] = 0.5
    sigma_th[('U233', 'Fission')] = 500.0
    sigma_fa[('U233', 'Fission')] = 1.5

    # --- Uranium 235 ---
    sigma_th[('U235', 'Fission')] = 580.0
    sigma_fa[('U235', 'Fission')] = 1.0
    sigma_th[('U235', 'Capture')] = 100.0
    sigma_fa[('U235', 'Capture')] = 0.3

    # --- U236 / U237 : essentiellement capture (résultats intermédiaires) ---
    sigma_th[('U236', 'Fission')] = 0.0
    sigma_fa[('U236', 'Fission')] = 0.1
    sigma_th[('U236', 'Capture')] = 5.0
    sigma_fa[('U236', 'Capture')] = 0.2

    sigma_th[('U237', 'Fission')] = 0.0
    sigma_fa[('U237', 'Fission')] = 0.1
    sigma_th[('U237', 'Capture')] = 2.0
    sigma_fa[('U237', 'Capture')] = 0.1

    # --- Uranium 238 ---
    sigma_th[('U238', 'Fission')] = 0.0      # quasi impossible thermique
    sigma_fa[('U238', 'Fission')] = 0.3
    sigma_th[('U238', 'Capture')] = 2.7
    sigma_fa[('U238', 'Capture')] = 0.3

    # --- U239 (intermédiaire de la chaîne fertile) ---
    sigma_th[('U239', 'Fission')] = 0.0
    sigma_fa[('U239', 'Fission')] = 0.1
    sigma_th[('U239', 'Capture')] = 5.0
    sigma_fa[('U239', 'Capture')] = 0.2

    # --- Neptunium 239 ---
    sigma_th[('Np239', 'Fission')] = 5.0
    sigma_fa[('Np239', 'Fission')] = 0.5
    sigma_th[('Np239', 'Capture')] = 50.0
    sigma_fa[('Np239', 'Capture')] = 0.5

    # --- Plutonium 239 ---
    sigma_th[('Pu239', 'Fission')] = 750.0
    sigma_fa[('Pu239', 'Fission')] = 1.8
    sigma_th[('Pu239', 'Capture')] = 270.0
    sigma_fa[('Pu239', 'Capture')] = 0.5

    # --- Plutonium 240 ---
    sigma_th[('Pu240', 'Fission')] = 0.0
    sigma_fa[('Pu240', 'Fission')] = 0.1
    sigma_th[('Pu240', 'Capture')] = 290.0
    sigma_fa[('Pu240', 'Capture')] = 0.5

    # --- Xénon 135 : poison neutronique très fort ---
    sigma_th[('Xe135', 'Fission')] = 0.0
    sigma_fa[('Xe135', 'Fission')] = 0.0
    sigma_th[('Xe135', 'Capture')] = 2.0e6   # énorme à 0.025 eV
    sigma_fa[('Xe135', 'Capture')] = 10.0    # beaucoup plus faible en rapide

    key = (X, Transfo)
    if key not in sigma_th:
        raise ValueError(f"No cross-section data for ({X}, {Transfo}) in this simplified model.")

    # ============================== Calcul sigma(E) ==============================

    sigma_thermal = sigma_th[key]
    sigma_fast = sigma_fa[key]

    # masque thermique / rapide
    thermal_mask = E < 0.5  # eV
    sigma = np.where(thermal_mask, sigma_thermal, sigma_fast)

    # ============================== Retour du bon type ==============================

    if is_scalar:
        return float(sigma[0])
    elif is_list:
        return sigma.tolist()
    else:
        return sigma


# Exemple de test manuel :
if __name__ == "__main__":
    energies = np.logspace(-5, 6, 5)
    cs = crossSection(X='U235', Transfo='Fission', E_neutron=energies)
    print("E (eV) :", energies)
    print("sigma_fission_U235 (barn) :", cs)
