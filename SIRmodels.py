__author__ = "Patricio Andrés Olivares Roncagliolo"
__email__ = "patricio.olivaresr@usm.cl"

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi

"""
    SIR class
"""
class SIR:
    def __init__(self, SIR0, params, t_sim):
        self.params = params
        self.SIR0 = SIR0
        self.t_sim = t_sim
        self.SIR_Res = None

    def __SIR_eqs(self, SIR0, t, params):
        # Initial conditions
        Si = SIR0[0]
        Ii = SIR0[1]
        Ri = SIR0[2]

        # Number of infected people
        N = np.sum(SIR0)

        # Parameters of SIR model
        beta, delta = params

        # Equations definition of the model
        SIR_S = - (beta * Si * Ii) / N
        SIR_I = (beta * Si * Ii) / N - delta * Ii
        SIR_R = delta * Ii

        return SIR_S, SIR_I, SIR_R

    def __modelSIR(self, SIR0, t, params):
        SIR_Res = spi.odeint(self.__SIR_eqs, SIR0, t, args=(params,))
        return SIR_Res

    def runSim(self):
        self.SIR_Res = self.__modelSIR(self.SIR0, self.t_sim, self.params)

    def getResult(self):
        return self.SIR_Res

    def plotSim(self):
        plt.plot(self.t_sim, self.SIR_Res[:, 1], '-r')
        plt.show()


class SIIR:
    def __init__(self, SIIR0, params, t_sim):
        self.SIIR0 = SIIR0
        self.params = params
        self.t_sim = t_sim
        self.SIIR_Res = None
        self.dis1 = None
        self.dis2 = None

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

        SIIR_SS = -SSi * beta1 * (ISi + IIi + IRi) / N - SSi * beta2 * (SIi + IIi + RIi) / N
        SIIR_IS = SSi * beta1 * (ISi + IIi + IRi) / N - delta1 * ISi - ISi * beta2prime * (SIi + IIi + RIi) / N
        SIIR_SI = SSi * beta2 * (SIi + IIi + RIi) / N - delta2 * SIi - SIi * beta1prime * (ISi + IIi + IRi) / N
        SIIR_II = ISi * beta2prime * (SIi + IIi + RIi) / N + SIi * beta1prime * (ISi + IIi + IRi) / N - delta1prime * IIi - delta2prime * IIi
        SIIR_RS = delta1 * ISi - RSi * beta2 * (SIi + IIi + RIi) / N
        SIIR_SR = delta2 * SIi - SRi * beta1 * (ISi + IIi + IRi) / N
        SIIR_RI = RSi * beta2 * (SIi + IIi + RIi) / N + delta1prime * IIi - delta2 * RIi
        SIIR_IR = SRi * beta1 * (ISi + IIi + IRi) / N + delta2prime * IIi - delta1 * IRi
        SIIR_RR = delta1 * IRi + delta2 * RIi

        return SIIR_SS, SIIR_IS, SIIR_SI, SIIR_II, SIIR_RS, SIIR_SR, SIIR_RI, SIIR_IR, SIIR_RR

    def __modelSIIR(self, SIIR0, t, params):
        SIIR_Res = spi.odeint(self.__SIIR_eqs, SIIR0, t, args=(params,))
        return SIIR_Res

    def runSim(self):
        self.SIIR_Res = self.__modelSIIR(self.SIIR0, self.t_sim, self.params)
        self.dis1 = self.SIIR_Res[:, 1] + self.SIIR_Res[:, 3] + self.SIIR_Res[:, 7]
        self.dis2 = self.SIIR_Res[:, 2] + self.SIIR_Res[:, 3] + self.SIIR_Res[:, 6]

    def getResult(self):
        return self.SIIR_Res

    def getDisease1(self):
        return self.dis1

    def getDisease2(self):
        return self.dis2

    def plotSim(self):
        fig, ax = plt.subplots(1, 2)
        ax[0].plot(self.t_sim, self.dis1, '-r')
        ax[1].plot(self.t_sim, self.dis2, '-b')
        plt.show()


def testSIR():
    beta = 14.3553504
    delta = 13.16275817
    params = (beta, delta)

    N = 21739040
    I0 = 119
    S0 = N - I0
    R0 = 0
    SIR0 = (S0,I0,R0)


    t_start = 0
    t_end = 12
    n_int = 10000

    t_sim = np.linspace(t_start, t_end, n_int)


    #SIR_Res, t = mdl.modelSIR(SIR0, t_sim, params)
    #plt.plot(t_sim, SIR_Res[:, 1], '-r')
    #plt.show()

    sirSim = SIR(SIR0, params, t_sim)
    sirSim.runSim()
    sirSim.plotSim()
    #print(sirSim.getResult())

def testSIIR():
    beta1 = 14.3553504
    beta1prime = 14.51848543
    delta1 = 13.16275817
    delta1prime = 13.98388944

    beta2 = 206.18920775
    beta2prime = 222.29093379
    delta2 = 206.04507902
    delta2prime = 154.4459061

    N = 21739040

    t_start = 0
    t_end = 12
    n_int = 10000

    t_sim = np.linspace(t_start, t_end, n_int)
    params = beta1, beta2, delta1, delta2, beta1prime, beta2prime, delta1prime, delta2prime

    SIIR0 = np.zeros(9)
    SIIR0[1] = 119
    SIIR0[2] = 7
    SIIR0[0] = N - np.sum(SIIR0[1:8])

    siirSim = SIIR(SIIR0, params, t_sim)
    siirSim.runSim()
    siirSim.plotSim()
    #print(siirSim.getResult())

