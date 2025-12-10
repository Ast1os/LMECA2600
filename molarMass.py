# molarMass.py
# Renvoie la masse molaire en kg/mol pour quelques noyaux utiles du projet.

def molarMass(nuclide: str) -> float:
    """
    Retourne la masse molaire (en kg/mol) du nucléide donné.

    Paramètres
    ----------
    nuclide : str
        Nom du nucléide sous forme "U235", "U238", "Pu239", "Xe135", etc.
        On peut aussi utiliser "n" pour le neutron.

    Retour
    ------
    float
        Masse molaire en kg/mol.

    Remarque
    --------
    Les masses sont approximées par A g/mol -> A/1000 kg/mol,
    ce qui est largement suffisant pour ce projet.
    """

    # masses molaires approximatives en g/mol
    molar_mass_g_per_mol = {
        # Thorium / Protactinium / Uranium (pour la chaîne fertile)
        "Th232": 232.0,
        "Th233": 233.0,
        "Pa233": 233.0,

        # Uranium
        "U235": 235.0,
        "U236": 236.0,
        "U238": 238.0,
        "U239": 239.0,

        # Neptunium
        "Np239": 239.0,

        # Plutonium
        "Pu239": 239.0,
        "Pu240": 240.0,

        # Xenon
        "Xe135": 135.0,

        # Neutron (masse ≈ 1 u)
        "n": 1.0,
    }

    if nuclide not in molar_mass_g_per_mol:
        raise ValueError(f"Molar mass not defined for nuclide '{nuclide}'")

    return molar_mass_g_per_mol[nuclide] / 1000.0  # conversion g/mol -> kg/mol


if __name__ == "__main__":
    # Petits tests rapides
    for X in ["U235", "U238", "Pu239", "Xe135", "n"]:
        print(X, molarMass(X), "kg/mol")
