__author__ = "Patricio Andrés Olivares Roncagliolo"
__email__ = "patricio.olivaresr@usm.cl"

import unittest

import numpy as np
from models.SIR import SIR
from models.exceptions import NotEvaluatedError


class TestSIRModel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Input parameters
        beta, gamma = 10, 5  # R0=2
        params = (beta, gamma)

        # Initial conditions
        N = 10000
        I0 = 100  # 0.1% infected people
        S0 = N - I0
        R0 = 0
        SIR0 = (S0, I0, R0)

        t_start = 0
        t_end = 100
        n_int = 10000

        t_sim = np.linspace(t_start, t_end, n_int)

        cls.sirSim = SIR(SIR0, params, t_sim)

    def test_NotEvaluatedError(self):
        # Tested values
        self.assertRaises(NotEvaluatedError, self.sirSim.getDisease)
        self.assertRaises(NotEvaluatedError, self.sirSim.getNInfected)
        self.assertRaises(NotEvaluatedError, self.sirSim.getPeakPos)
        self.assertRaises(NotEvaluatedError, self.sirSim.getResult)


def testSIR():
    beta = 5.06
    delta = 5
    params = (beta, delta)

    N = 1000000
    I0 = 100
    S0 = N - I0
    R0 = 0
    SIR0 = (S0, I0, R0)

    t_start = 0
    t_end = 100
    n_int = 10000

    t_sim = np.linspace(t_start, t_end, n_int)

    # SIR_Res, t = mdl.modelSIR(SIR0, t_sim, params)
    # plt.plot(t_sim, SIR_Res[:, 1], '-r')
    # plt.show()

    sirSim = SIR(SIR0, params, t_sim)
    sirSim.runEvaluation()
    res = sirSim.getResult()
    plt.plot(t_sim, res[:, 0])
    plt.plot(t_sim, res[:, 1])
    plt.plot(t_sim, res[:, 2])
    plt.show()
    print(sirSim.getDisease()[1])
    # sirSim.plotSeries()
    # print(t_sim[sirSim.getResult()[1]])
    # print(sirSim.getResult())
