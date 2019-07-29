__author__ = "Patricio Andrés Olivares Roncagliolo"
__email__ = "patricio.olivaresr@usm.cl"

import SIRmodels as mdl
import numpy as np
import matplotlib.pyplot as plt

class modelAnalysisSIR:
    def __init__(self,SIR0, params, t_sim):
        self.sir = mdl.SIR(SIR0, params, t_sim)
        self.params = params
        self.t_sim = t_sim
        self.SIR0 = SIR0


    def plotBetaPeak(self, beta_arr=None, savefig=False, norm=False):
        if beta_arr is None:
            beta_arr = np.array([self.params[0]])

        params = self.params.copy()

        peak = np.zeros(len(beta_arr))
        for i in np.arange(len(beta_arr)):
            params[0] = beta_arr[i]
            sir = mdl.SIR(self.SIR0, params, self.t_sim)
            sir.runEvaluation(norm)
            try:
                peak[i] = self.t_sim[sir.getDisease()[1]][0]
            except:
                peak[i] = self.t_sim[0]

        plt.plot(beta_arr, peak)
        plt.xlabel("Beta")
        plt.ylabel("Peak position")
        txt = ("Beta = " + str(self.params[0]) + "\nDelta = " + str(self.params[1]))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        plt.text(0.6, 0.4, txt, bbox=props)
        plt.show()




def testAnalysis():
    # Originals values
    N = 1000000
    SIR0 = np.zeros(3)
    SIR0[1] = 100
    SIR0[0] = N - np.sum(SIR0[1:2])

    t_start = 0
    t_end = 100
    n_int = 10000

    t_sim = np.linspace(t_start, t_end, n_int)
    params = [10, 5]

    sirSim_orig = modelAnalysisSIR(SIR0, params, t_sim)
    sirSim_orig.plotBetaPeak(beta_arr=np.arange(0, 40, 0.01))

def getAnalysis():
    # Originals values
    N = 1000000
    SIIR0 = np.zeros(9)
    SIIR0[1] = 100
    SIIR0[2] = 100
    SIIR0[0] = N - np.sum(SIIR0[1:8])

    t_start = 0
    t_end = 15
    n_int = 10000

    t_sim = np.linspace(t_start, t_end, n_int)
    #params = [10, 10, 5, 9.9, 10, 10, 3, 8]
    #params = [10, 10, 5, 10, 10, 10, 5, 5]
    #params = [20, 10, 15, 9, 20, 10, 19.5, 9]

    params = [10, 10, 9, 5, 10, 10, 2, 5] # Fixed delta primes
    #params = [10, 10, 5, 9, 10, 10, 5, 5]  # Fixed delta primes


    siirSim_orig = modelAnalysis(SIIR0, params, t_sim)

    siirSim_orig.plotChange(betaprime_arr=np.array([5, 15, 40]), fixed=1, peak=False, norm=True)
    siirSim_orig.plotChange(betaprime_arr=np.array([5, 15, 40]), fixed=2, peak=False, norm=True)
    siirSim_orig.plotPeak2D(beta1prime_arr=np.arange(0, 40, 2), beta2prime_arr=np.arange(0, 40, 2), savefig=False, norm=True)
    siirSim_orig.plotInfected2D(beta1prime_arr=np.arange(0, 40, 2), beta2prime_arr=np.arange(0, 40, 2), savefig=False, norm=True)
    siirSim_orig.plotPeakNInfected(beta1prime_arr=np.arange(0, 40, 2), beta2prime_arr=np.arange(0, 40, 2), savefig=False, norm=True)

testAnalysis()