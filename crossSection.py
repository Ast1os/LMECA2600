import numpy as np


def crossSection(X, Transfo, E_neutron):
    """
    Computes the cross section for a specific energy level
    ENDF database: https://www-nds.iaea.org/exfor/endf.htm

    ----------------
    :param X: string
        nuclide or nucleon (for a neutron), that follows the atomic notation of the element,
        e.g.: for the Uranium 235, X = 'U235'
    :param Transfo: string
        name of the considered transformation, should be one of the following in the list
        ['Fission', 'Capture']
    :param E_neutron: array, list
        value(s) of the cross section for an incident neutron of energy level equal to E_neutron in [eV]
        the value should be in the range [1e-5; 2e7] [eV]
    :return: array, list
        values of the cross section in [barn] corresponding to the input parameters
    """

    if X != 'Th232' and X != 'Th233' and X != 'Pa233' and X != 'U233' and X != 'U235' and X != 'U236' and X != 'U237' and X != 'U238' \
            and X != 'U239' and X != 'Np239' and X != 'Pu239' and X != 'Pu240' and X != 'Xe135':
        print('\n WARNING : There is no database for element ', X, '. \n Please check function information')

    E = np.array(E_neutron, dtype=float)

    if np.any(E <= 0.0):
        raise ValueError("E_neutron must be > 0 eV")

    cs_data = {
        'Th232': {
            'Capture': {'thermal': 7.4,  'fast': 0.3},
            'Fission': {'thermal': 0.0,  'fast': 0.2},
        },
        'Th233': {
            'Capture': {'thermal': 20.0, 'fast': 0.5},
            'Fission': {'thermal': 0.0,  'fast': 0.1},
        },
        'Pa233': {
            'Capture': {'thermal': 140.0, 'fast': 5.0},
            'Fission': {'thermal': 0.0,   'fast': 0.1},
        },
        'U233': {
            'Capture': {'thermal': 45.0,  'fast': 0.05},
            'Fission': {'thermal': 531.0, 'fast': 1.0},
        },
        'U235': {
            'Capture': {'thermal': 99.0,  'fast': 0.09},
            'Fission': {'thermal': 583.0, 'fast': 1.0},
        },
        'U236': {
            'Capture': {'thermal': 6.0,   'fast': 0.1},
            'Fission': {'thermal': 0.1,   'fast': 0.1},
        },
        'U237': {
            'Capture': {'thermal': 10.0,  'fast': 0.2},
            'Fission': {'thermal': 0.0,   'fast': 0.1},
        },
        'U238': {
            'Capture': {'thermal': 2.0,   'fast': 0.07},
            'Fission': {'thermal': 2e-5,  'fast': 0.3},
        },
        'U239': {
            'Capture': {'thermal': 10.0,  'fast': 0.2},
            'Fission': {'thermal': 0.0,   'fast': 0.1},
        },
        'Np239': {
            'Capture': {'thermal': 20.0,  'fast': 0.5},
            'Fission': {'thermal': 0.5,   'fast': 0.5},
        },
        'Pu239': {
            'Capture': {'thermal': 269.0, 'fast': 0.05},
            'Fission': {'thermal': 748.0, 'fast': 2.0},
        },
        'Pu240': {
            'Capture': {'thermal': 290.0, 'fast': 0.1},
            'Fission': {'thermal': 0.05,  'fast': 0.2},
        },
        'Xe135': {
            'Capture': {'thermal': 2.0e6, 'fast': 8.0e-4},
            'Fission': {'thermal': 0.0,   'fast': 0.0},
        },
    }

    if X not in cs_data:
        E_shape = E.shape
        sigma = np.zeros_like(E)
        if sigma.shape == () and not np.iterable(E_neutron):
            return float(sigma)
        return sigma

    if Transfo not in cs_data[X]:
        E_shape = E.shape
        sigma = np.zeros_like(E)
        if sigma.shape == () and not np.iterable(E_neutron):
            return float(sigma)
        return sigma

    E_boundary = 1.0

    sigma_th = cs_data[X][Transfo]['thermal']
    sigma_fast = cs_data[X][Transfo]['fast']

    sigma = np.where(E < E_boundary, sigma_th, sigma_fast)

    if sigma.shape == () and not np.iterable(E_neutron):
        return float(sigma)

    return sigma