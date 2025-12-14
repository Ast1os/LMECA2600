import numpy as np


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

    half_life_db = {
        'Th232': {
            'Alpha': 4.4338428e17,
        },
        'Th233': {
            'BetaMinus': 1338,
        },
        'Pa233': {
            'BetaMinus': 2.33064e6,
        },
        'U233': {
            'Alpha': 5.02396992e12,
        },
        'U235': {
            'Alpha': 2.221950368e16,
        },
        'U236': {
            'Alpha': 7.39078992e14,
        },
        'U237': {
            'BetaMinus': 5.832e5,
        },
        'U238': {
            'Alpha': 1.409993568e17,
        },
        'U239': {
            'BetaMinus': 1407.0,
        },
        'Np239': {
            'BetaMinus': 2.035584e5,
        },
        'Pu239': {
            'Alpha': 7.60853736e11,
        },
        'Pu240': {
            'Alpha': 2.070494136e11,
        },
        'Xe135': {
            'BetaMinus': 3.2904e4,
        },
    }

    if X != 'Th232' and X != 'Th233' and X != 'Pa233' and X != 'U233' and X != 'U235'and X != 'U236'and X != 'U237' and X != 'U238'\
            and X != 'U239' and X != 'Np239' and X != 'Pu239' and X != 'Pu240' and X != 'Xe135':
        print('\n WARNING : There is no database for element ', X, '. \n Please check function information')

    if Transfo != 'Alpha' and Transfo != 'BetaMinus' and Transfo != 'BetaPlus' and Transfo != 'Gamma':
        print('\n WARNING : These transformation are not implemented :', Transfo, '.\n Please check function information')

    if X in half_life_db and Transfo in half_life_db[X]:
        hl = half_life_db[X][Transfo]
    else:
        print('\n WARNING : No half-life data for', X, 'with transformation', Transfo,
              '. Assuming effectively stable (hl = np.inf).\n')
        hl = np.inf

    return hl
