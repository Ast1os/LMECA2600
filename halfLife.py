# halfLife.py
# Renvoie la demi-vie (en secondes) pour quelques noyaux et modes de décroissance.

def halfLife(nuclide: str, transfo: str) -> float:
    """
    Retourne la demi-vie T1/2 (en secondes) du nucléide 'nuclide'
    pour un mode de transformation radioactif donné.

    Paramètres
    ----------
    nuclide : str
        Par ex. "U235", "U238", "U239", "Np239", "Pu239", "Pu240", "Xe135", etc.
    transfo : str
        Type de transformation radioactives, parmi :
        "Alpha", "BetaMinus", "BetaPlus", "Gamma"

    Retour
    ------
    float
        Demi-vie en secondes.

    Exceptions
    ----------
    ValueError si la combinaison (nuclide, transfo) n'est pas définie.

    Remarque
    --------
    Les valeurs numériques sont approximatives mais physiquement réalistes.
    Pour ce projet, seules les demi-vies à l'échelle de la seconde / heure
    jouent un rôle (Xe135, précurseurs...), les demi-vies géologiques peuvent
    être considérées comme infinies à l'échelle du transitoire.
    """

    # Constantes pratiques
    year = 365.25 * 24 * 3600.0
    hour = 3600.0
    minute = 60.0
    day = 24 * 3600.0

    # dictionnaire : (nuclide, transfo) -> T1/2 en s
    hl = {
        # Chaîne fertile U238 -> Pu239
        ("U239", "BetaMinus"): 23.5 * minute,     # ~23.5 min
        ("Np239", "BetaMinus"): 2.356 * day,      # ~2.36 jours

        # Xe135 (important pour l'effet poison)
        ("Xe135", "BetaMinus"): 9.14 * hour,      # ~9.14 h

        # Actinides lourds : demi-vies énormes (quasi stables pour nous)
        ("U235", "Alpha"): 7.04e8 * year,
        ("U238", "Alpha"): 4.47e9 * year,
        ("Th232", "Alpha"): 1.40e10 * year,
        ("Pu239", "Alpha"): 2.41e4 * year,
        ("Pu240", "Alpha"): 6.56e3 * year,

        # Chaîne fertile Th232 -> U233 (au cas où)
        ("Th233", "BetaMinus"): 22.3 * minute,
        ("Pa233", "BetaMinus"): 26.97 * day,
    }

    key = (nuclide, transfo)
    if key not in hl:
        raise ValueError(f"Half-life not defined for nuclide '{nuclide}' with transformation '{transfo}'")

    return hl[key]


if __name__ == "__main__":
    tests = [
        ("U239", "BetaMinus"),
        ("Np239", "BetaMinus"),
        ("Xe135", "BetaMinus"),
        ("U235", "Alpha"),
    ]
    for X, T in tests:
        print(X, T, ":", halfLife(X, T), "s")

