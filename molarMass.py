def molarMass(X):
    """
    Database: http://wwwndc.jaea.go.jp/NuC/
    :param X: string
        nuclide or nucleon (for a neutron), that follows the atomic notation of the element,
        i.e.: for the Uranium 235, X = 'U235'; for a neutron X = 'n'
    :return: double
        molar mass of the nucleon, or nuclide X in [kg/mol]
    """

    # Your work : add all species

    if X == 'Th232':
        M = 232.038060026e-3
    elif X == 'Th233':
        M = 233.041586541e-3
    elif X == 'Pa233':
        M = 233.040248815e-3
    elif X == 'U233':
        M = 233.039636574e-3
    elif X == 'U235':
        M = 235.043931368e-3
    elif X == 'U236':
        M = 236.045569468e-3
    elif X == 'U237':
        M = 237.048731636e-3
    elif X == 'U238':
        M = 238.050789466e-3
    elif X == 'U239':
        M = 239.054294518e-3
    elif X == 'Np239':
        M = 239.052940487e-3

    elif X == 'Pu239':
        M = 239.052164844e-3
    elif X == 'Pu240':
        M = 240.053815008e-3

    elif X == 'Xe135':
        M = 134.907226844e-3
    # Neutron
    elif X == 'n':
        # Neutron molar mass (CODATA 2018): ~1.008664917×10^-3 kg/mol
        M = 1.008664917e-3

    else:
        print('\n ---------------- WARNING -----------------\n No molar mass for species (', X, ')')
        M = 0.0

    return M

# Try your function:
Mm = molarMass('U235')
print(Mm)