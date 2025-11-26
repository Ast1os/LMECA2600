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
        M = 0
    elif X == 'Th233':
        M = 0
    else:
        print('\n ---------------- WARNING ----------------- \n No molar mass for species (', X, ')')
        M = 0

    return M

# Try your function:
# Mm = molarMass('U235')
# print(Mm)