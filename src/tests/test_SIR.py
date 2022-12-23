__author__ = "Patricio Andrés Olivares Roncagliolo"
__email__ = "patricio.olivaresr@usm.cl"

import unittest

import numpy as np
from models.SIR import SIR
from models.exceptions import NotEvaluatedError


class TestSIRModelMethods(unittest.TestCase):
    def test_case1(self):
        # Input parameters
        beta, gamma = 2, 0
        params = (beta, gamma)

        # Initial conditions
        N = 10000
        I0 = 100  # 1% infected people
        S0 = N - I0
        R0 = 0
        SIR0 = (S0, I0, R0)

        t_start = 0
        t_end = 1000
        n_int = 10000

        t_sim = np.linspace(t_start, t_end, n_int)

        sirSim = SIR(SIR0, params, t_sim)
        sirSim.run()

        # Testing infected number
        self.assertEqual(10000, sirSim.getNInfected())

        SIR_res = sirSim.getResult()
        # Random picking 100 values for evaluating if sumation equal 10000 (entire population)
        selection_index = np.random.randint(0, t_sim.size, 100)
        for i in selection_index:
            count_sum = SIR_res[i, 0] + SIR_res[i, 1] + SIR_res[i, 2]
            self.assertAlmostEqual(10000, count_sum, delta=2)

        # Should not be peak position
        peakpos = sirSim.getPeakPos()
        self.assertTrue(len(peakpos) == 0)

    def test_case2(self):
        # Input parameters
        beta, gamma = 2, 1  # R0 2
        params = (beta, gamma)

        # Initial conditions
        N = 10000
        I0 = 100  # 1% infected people
        S0 = N - I0
        R0 = 0
        SIR0 = (S0, I0, R0)

        t_start = 0
        t_end = 1000
        n_int = 10000

        t_sim = np.linspace(t_start, t_end, n_int)

        sirSim = SIR(SIR0, params, t_sim)
        sirSim.run()

        # Testing infected number
        self.assertGreaterEqual(10000, sirSim.getNInfected())

        SIR_res = sirSim.getResult()
        # Random picking 100 values for evaluating if sumation equal 10000 (entire population)
        selection_index = np.random.randint(0, t_sim.size, 100)
        for i in selection_index:
            count_sum = SIR_res[i, 0] + SIR_res[i, 1] + SIR_res[i, 2]
            self.assertAlmostEqual(10000, count_sum, delta=2)

    def test_case3(self):
        # Input parameters
        beta, gamma = 1, 2  # R0 0.5
        params = (beta, gamma)

        # Initial conditions
        N = 10000
        I0 = 100  # 1% infected people
        S0 = N - I0
        R0 = 0
        SIR0 = (S0, I0, R0)

        t_start = 0
        t_end = 1000
        n_int = 10000

        t_sim = np.linspace(t_start, t_end, n_int)

        sirSim = SIR(SIR0, params, t_sim)
        sirSim.run()

        # Testing infected number
        self.assertGreaterEqual(10000, sirSim.getNInfected())

        SIR_res = sirSim.getResult()
        # Random picking 100 values for evaluating if sumation equal 10000 (entire population)
        selection_index = np.random.randint(0, t_sim.size, 100)
        for i in selection_index:
            count_sum = SIR_res[i, 0] + SIR_res[i, 1] + SIR_res[i, 2]
            self.assertAlmostEqual(10000, count_sum, delta=2)
            self.assertGreaterEqual(100, SIR_res[i, 1])

        # Should not be peak position
        peakpos = sirSim.getPeakPos()
        self.assertTrue(len(peakpos) == 0)

    def test_NotEvaluatedError(self):
        # Input parameters
        beta, gamma = 10, 5  # R0=2
        params = (beta, gamma)

        # Initial conditions
        N = 9900
        I0 = 100  # 1% infected people
        S0 = N - I0
        R0 = 0
        SIR0 = (S0, I0, R0)

        t_start = 0
        t_end = 100
        n_int = 10000

        t_sim = np.linspace(t_start, t_end, n_int)

        sirSim = SIR(SIR0, params, t_sim)
        # Tested values
        self.assertRaises(NotEvaluatedError, sirSim.getDisease)
        self.assertRaises(NotEvaluatedError, sirSim.getNInfected)
        self.assertRaises(NotEvaluatedError, sirSim.getPeakPos)
        self.assertRaises(NotEvaluatedError, sirSim.getResult)