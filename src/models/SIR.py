__author__ = "Patricio Andrés Olivares Roncagliolo"
__email__ = "patricio.olivaresr@usm.cl"

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import peakutils as pk


class SIR:
    """ SIR Class that run the original SIR model

    Attributes:
        SIR0 (list): Initial conditions for the SIR model.
        params (list): Parameters for the SIR model. It includes the infection rate and recovery rate.
        t_sim (list): List that includes the time range of the simulation.
        SIR_Res (list): Saves the result of the running SIR model (time series) for the given conditions for all the states.
        peakpos (list): Include the list position of all peaks of the time series result.

    """

    def __init__(self, SIR0, params, t_sim) -> None:
        """ Constructor

        :param SIR0: Initial conditions for the SIR model.
        :type SIR0: list
        :param params: Parameters for the SIR model. It includes the infection rate and recovery rate.
        :type params: list
        :param t_sim: List that includes the time range of the simulation.
        :type t_sim: list
        """
        self.SIR0 = SIR0
        self.params = params
        self.t_sim = t_sim
        self.SIR_res = None
        self.peakpos = None

    def __SIR_eqs(self, SIR0, t, params):
        """ Private method that include the SIR model equations

        :param SIR0: Initial conditions for the SIR model.
        :type SIR0: list
        :param t: List that includes the time range of the simulation.
        :type t: list
        :param params: Parameters for the SIR model. It includes the infection rate and recovery rate.
        :type params: list
        :return: The result for the Susceptible (S), Infected (I) and Recovery (R) states
        :rtype: list
        """

        # Initial conditions
        Si = SIR0[0]
        Ii = SIR0[1]
        Ri = SIR0[2]

        # Number of infected people
        N = np.sum(SIR0)

        # Parameters of SIR model
        beta, delta = params

        # Equations definition of the model
        S = - (beta * Si * Ii) / N
        I = (beta * Si * Ii) / N - delta * Ii
        R = delta * Ii

        return S, I, R

    # Configure SIR model with initial conditions, parameters and time
    def __modelSIR(self, SIR0, t, params):
        """ Private method that configure the SIR model with initial conditions, parameters and time.

        :param SIR0: Initial conditions for the SIR model.
        :type SIR0: list
        :param t: List that includes the time range of the simulation.
        :type t: list
        :param params: Parameters for the SIR model. It includes the infection rate and recovery rate.
        :type params: list
        :return: Final result of the evaluation using the SIR equations (SIR_res)
        :rtype: 2D-list
        """
        SIR_res = spi.odeint(self.__SIR_eqs, SIR0, t, args=(params,))
        return SIR_res

    def run(self, norm=False):
        """
        Run the evaluation of the model and saves the result on SIR_Res variable and the peak position of the infection
        of the disease on peakpos variable.
        :return: SIR_Res
        """
        self.SIR_res = self.__modelSIR(self.SIR0, self.t_sim, self.params)
        if norm:
            N = np.sum(self.SIR0)
            self.SIR_res = self.SIR_res / N
            
        # Check for peak infection position and save it
        self.peakpos = pk.indexes(self.SIR_tes[:, 1], thres=0.5)

    def getResult(self):
        """
        Return the result of the evaluation (SIR_Res) and the peak position of the infection (peakpos)
        :return: SIR_Res
        """
        return self.SIR_res

    def getDisease(self):
        return self.SIR_res[:, 1], self.peakpos

    def getNInfected(self):
        n_infected = self.SIR_res[:, 2][-1]
        return n_infected

    # Plot time series result
    def plotSeries(self):
        """
        Plot the time series of the infected state over time.
        :return:
        """
        plt.plot(self.t_sim, self.SIR_res[:, 1], '-r')
        plt.show()




def testSIIR():

    beta1 = 5
    beta1prime = 5
    delta1 = 1
    delta1prime = 5

    beta2 = 10
    beta2prime = 20
    delta2 = 9.9
    delta2prime = 9.9

    N = 1000

    t_start = 0
    t_end = 10
    n_int = 10000

    t_sim = np.linspace(t_start, t_end, n_int)
    params = beta1, beta2, delta1, delta2, beta1prime, beta2prime, delta1prime, delta2prime

    SIIR0 = np.zeros(9)
    SIIR0[1] = 1
    SIIR0[2] = 1
    SIIR0[0] = N - np.sum(SIIR0[1:8])

    siirSim = SIIR(SIIR0, params, t_sim)
    siirSim.runEvaluation(norm=True)
    # siirSim.plotSeries()
    # siirSim.plotDisease1Series(savefig=True)
    # siirSim.plotDisease2Series(savefig=False)
    print("N infected1: " + str(siirSim.getNInfected1()))
    res = siirSim.getResult()
    # dy = np.trapz(res[:, 1] + res[:, 3] + res[:, 7], t_sim)
    # print("N infected1: " + str(np.sum(dy)))
    # print("N infected2: " + str(siirSim.getNInfected2()))
    plt.plot(t_sim, res[:, 0], 'g')
    plt.plot(t_sim, res[:, 1] + res[:, 3] + res[:, 6], 'r')
    plt.plot(t_sim, res[:, 4], 'b')
    plt.show()
