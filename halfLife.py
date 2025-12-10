def halfLife(X, Transfo):
    """
    Return the half life of a nuclide or nucleon for a specific transformation
    Data available on http://wwwndc.jaea.go.jp/NuC/ (more precisely: http://wwwndc.jaea.go.jp/CN14/index.html)

    ----------------
    :param X: string
        nuclide or nucleon (for a neutron), that follows the atomic notation of the element,
        i.e.: for the Uranium 235, X = 'U235'
    :param Transfo: string
        name of the considered transformation, should be one of the following in the list
        ['Alpha', 'BetaMinus', 'BetaPlus', 'Gamma']
    :return: the half life of X for the transformation Transfo in seconds
    """

    # =====================  Base de données des demi-vies  =====================
    # Toutes les valeurs sont en secondes.
    # Les très longues demi-vies (U-238, Pu-239, ...) sont gardées pour info,
    # mais seront en pratique quasi-stables sur l’échelle de temps du réacteur.

    half_life_db = {
        'Th232': {  # 1.405×10^10 years, alpha decay
            'Alpha': 4.4338428e17,   # s  :contentReference[oaicite:0]{index=0}
        },
        'Th233': {  # 21.83 min, beta-
            'BetaMinus': 1338,     # s  :contentReference[oaicite:1]{index=1}
        },
        'Pa233': {  # 26.975 days, beta-
            'BetaMinus': 2.33064e6,  # s  :contentReference[oaicite:2]{index=2}
        },
        'U233': {   # 1.592×10^5 years, alpha
            'Alpha': 5.02396992e12,  # s  :contentReference[oaicite:3]{index=3}
        },
        'U235': {   # 7.038×10^8 years, alpha
            'Alpha': 2.221950368e16,  # s  :contentReference[oaicite:4]{index=4}
        },
        'U236': {   # 2.342×10^7 years, alpha
            'Alpha': 7.39078992e14,  # s  :contentReference[oaicite:5]{index=5}
        },
        'U237': {   # ~6.75 days, beta-
            'BetaMinus': 5.832e5,    # s  :contentReference[oaicite:6]{index=6}
        },
        'U238': {   # 4.468×10^9 years, alpha
            'Alpha': 1.409993568e17, # s  :contentReference[oaicite:7]{index=7}
        },
        'U239': {   # 23.45 min, beta-
            'BetaMinus': 1407.0,     # s  :contentReference[oaicite:8]{index=8}
        },
        'Np239': {  # 2.356 days, beta-
            'BetaMinus': 2.035584e5, # s  :contentReference[oaicite:9]{index=9}
        },
        'Pu239': {  # 24110 years, alpha
            'Alpha': 7.60853736e11,  # s  :contentReference[oaicite:10]{index=10}
        },
        'Pu240': {  # 6561 years, alpha
            'Alpha': 2.070494136e11, # s  :contentReference[oaicite:11]{index=11}
        },
        'Xe135': {  # 9.14 h, beta-
            'BetaMinus': 3.2904e4,   # s  :contentReference[oaicite:12]{index=12}
        },
        # Tu peux éventuellement ajouter ici tes "FP" lumpés avec T1/2 ≈ 1 s,
        # par ex. 'FP': {'BetaMinus': 1.0}
    }

    # ==================================  Check arguments  ==================================
    #  Check that the nucleon/nuclide asked, and that the associated transformation exists in the database
    if X != 'Th232' and X != 'Th233' and X != 'Pa233' and X != 'U233' and X != 'U235'and X != 'U236'and X != 'U237' and X != 'U238'\
            and X != 'U239' and X != 'Np239' and X != 'Pu239' and X != 'Pu240' and X != 'Xe135':
        print('\n WARNING : There is no database for element ', X, '. \n Please check function information')

     # Check whether the transformation exists or not
    if Transfo != 'Alpha' and Transfo != 'BetaMinus' and Transfo != 'BetaPlus' and Transfo != 'Gamma': # Transfos Gamma très courtes
        print('\n WARNING : These transformation are not implemented :', Transfo, '.\n Please check function information')

    # =========================  Récupération de la demi-vie  =========================
    if X in half_life_db and Transfo in half_life_db[X]:
        hl = half_life_db[X][Transfo]
    else:
        # Si la transformation n’existe pas pour ce noyau, on le traite comme
        # "effectivement stable" pour le modèle (pas de décroissance sur les temps
        # considérés). Ça évite les divisions par zéro sur lambda = ln(2)/hl.
        print('\n WARNING : No half-life data for', X, 'with transformation', Transfo,
              '. Assuming effectively stable (hl = np.inf).\n')
        hl = np.inf

    return hl

# Try your function
# hl = halfLife(X='Pa233', Transfo='BetaMinus')
