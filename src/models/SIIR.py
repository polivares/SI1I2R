__author__ = "Patricio Andrés Olivares Roncagliolo"
__email__ = "patricio.olivaresr@usm.cl"

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import peakutils as pk


class SIIR:
    """ SIIR class
    """

    def __init__(self, SIIR0, params, t_sim):
        self.SIIR0 = SIIR0  # Initial conditions
        self.params = params  # SIIR configuration parameters: beta1, beta2, delta1, delta2, beta1prime, beta2prime, delta1prime, delta2prime
        self.t_sim = t_sim  # Time of simulation
        self.SIIR_Res = None  # Results of Simulation
        self.dis1 = None  # Time series result of disease one
        self.dis2 = None  # Time series result of disease two
        self.peak1pos = None
        self.peak2pos = None

    # SIIR equations
    def __SIIR_eqs(self, SIIR0, t, params):
        SSi = SIIR0[0]
        ISi = SIIR0[1]
        SIi = SIIR0[2]
        IIi = SIIR0[3]
        RSi = SIIR0[4]
        SRi = SIIR0[5]
        RIi = SIIR0[6]
        IRi = SIIR0[7]
        RRi = SIIR0[8]

        N = np.sum(SIIR0)

        beta1, beta2, delta1, delta2, beta1prime, beta2prime, delta1prime, delta2prime = params

        SS = -SSi * beta1 * (ISi + IIi + IRi) / N - SSi * \
            beta2 * (SIi + IIi + RIi) / N
        IS = SSi * beta1 * (ISi + IIi + IRi) / N - delta1 * \
            ISi - ISi * beta2prime * (SIi + IIi + RIi) / N
        SI = SSi * beta2 * (SIi + IIi + RIi) / N - delta2 * \
            SIi - SIi * beta1prime * (ISi + IIi + IRi) / N
        II = ISi * beta2prime * (SIi + IIi + RIi) / N + SIi * beta1prime * (
            ISi + IIi + IRi) / N - delta1prime * IIi - delta2prime * IIi
        RS = delta1 * ISi - RSi * beta2 * (SIi + IIi + RIi) / N
        SR = delta2 * SIi - SRi * beta1 * (ISi + IIi + IRi) / N
        RI = RSi * beta2 * (SIi + IIi + RIi) / N + \
            delta1prime * IIi - delta2 * RIi
        IR = SRi * beta1 * (ISi + IIi + IRi) / N + \
            delta2prime * IIi - delta1 * IRi
        RR = delta1 * IRi + delta2 * RIi

        return SS, IS, SI, II, RS, SR, RI, IR, RR

    # Configure SIIR model with initial conditions, parameters and time
    def __modelSIIR(self, SIIR0, t, params):
        SIIR_Res = spi.odeint(self.__SIIR_eqs, SIIR0, t, args=(params,))
        return SIIR_Res

    # Run simulation
    def runEvaluation(self, norm=False):
        self.SIIR_Res = self.__modelSIIR(self.SIIR0, self.t_sim, self.params)
        N = np.sum(self.SIIR0)
        if norm:
            self.SIIR_Res = self.SIIR_Res/N
        self.dis1 = self.SIIR_Res[:, 1] + \
            self.SIIR_Res[:, 3] + self.SIIR_Res[:, 7]
        self.dis2 = self.SIIR_Res[:, 2] + \
            self.SIIR_Res[:, 3] + self.SIIR_Res[:, 6]
        self.peak1pos = pk.indexes(self.dis1, thres=0.9)
        self.peak2pos = pk.indexes(self.dis2, thres=0.9)

    # Get Time Series of both diseases (whole result)
    def getResult(self):
        return self.SIIR_Res

    # Get Time Series of Disease one
    def getDisease1(self):
        return self.dis1, self.peak1pos

    # Get Time Series of Disease two
    def getDisease2(self):
        return self.dis2, self.peak2pos

    def getNInfected1(self):
        n_infected = self.SIIR_Res[:, 4][-1] + \
            self.SIIR_Res[:, 6][-1] + self.SIIR_Res[:, 8][-1]
        return n_infected

    def getNInfected2(self):
        n_infected = self.SIIR_Res[:, 5][-1] + \
            self.SIIR_Res[:, 7][-1] + self.SIIR_Res[:, 8][-1]
        return n_infected

    # Plot both Series
    def plotSeries(self):
        fig, ax = plt.subplots(1, 2)
        ax[0].plot(self.t_sim, self.dis1, '-r')
        ax[1].plot(self.t_sim, self.dis2, '-b')
        plt.show()

    def plotDisease1Series(self, peak=True, savefig=False):
        fig, ax = plt.subplots()
        ax.plot(self.t_sim, self.dis1, '-r')
        legends = ["Disease 1"]
        if peak:
            ax.plot(self.t_sim[self.peak1pos], self.dis1[self.peak1pos], '*b')
            legends.append("Peak(s) Position")
        ax.legend(legends)
        plt.suptitle("Disease 1")
        if savefig:
            fig.savefig("images/disease1.svg", format='svg')
        else:
            plt.show()

    def plotDisease2Series(self, peak=True, savefig=False):
        fig, ax = plt.subplots()
        ax.plot(self.t_sim, self.dis2, '-r')
        legends = ["Disease 2"]
        if peak:
            ax.plot(self.t_sim[self.peak2pos], self.dis2[self.peak2pos], '*b')
            legends.append("Peak(s) Position")
        ax.legend(legends)
        plt.suptitle("Disease 2")
        if savefig:
            fig.savefig("images/disease2.svg", format='svg')
        else:
            plt.show()


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
