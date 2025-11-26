import numpy as np
import sys
import matplotlib.pyplot as plt

import molarMass as mM
import halfLife as hL
import crossSection as cS


def reactorModel(fuelCompo, FPCompo, t_final, n_th_init, n_fa_init, mTot):

    """
    Computes the evolution of the species in the reactor, and compute the power evolution in time
    :param fuelCompo: class
    :param FPCompo: class
            FPCompo.Xe135: percentage of FP that are considered as poison (in this project, we assume only
                Xe135 is produced as a poison)
            FPCompo.FP: percentage of FP that are considered "others fission products" of U235 (assume only one isotope
                different than Xe135 is produced by the fission)
    :param t_final: double
            final time of the simulation in seconds [s]
    :param n_th_init: double
            initial number of thermal neutrons
    :param n_fa_init: double
            initial number of fast neutrons
    :param mTot: double
            total number of kg of fuel at the initial state in [kg]
    :return fuelBurnup: array
            depending on the fuel, returns the burnup
    """

    return fuelBurnup

class Fuel:
    def __init__(self):
        self.U235 = 3
        self.U238 = 97
        self.Pu239 = 0
        self.Th232 = 0

class FP:
    def __init__(self):
        self.Xe135 = 5
        self.FP = 95

fuelCompo = Fuel()
FPCompo = FP()


# Try the function

# Mm = mM.molarMass('U235')
# hl = hL.halfLife(X='Pa233', Transfo='BetaMinus')
# cs = cS.crossSection(X='Pu240', Transfo='Fission', E_neutron=np.logspace(-5,6,10000).tolist())
#
# print(Mm, hl, cs)

# fBurnup = reactorModel(fuelCompo=fuelCompo, FPCompo=FPCompo, t_final=200., n_th_init=1e+10, n_fa_init=0., mTot=25.)