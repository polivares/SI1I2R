__author__ = "Patricio Andrés Olivares Roncagliolo"
__email__ = "patricio.olivaresr@usm.cl"

import SIRmodels as mdl
import numpy as np
import matplotlib.pyplot as plt

class modelAnalysis:
    def __init__(self,SIIR0, params, t_sim):
        self.siir = mdl.SIIR(SIIR0, params, t_sim)
        self.params = params
        self.t_sim = t_sim
        self.SIIR0 = SIIR0

    #
    # fixed: Disease fixed
    def plotParamsChange(self, beta1prime_arr=None, beta2prime_arr=None, fixed=1):
        if beta1prime_arr is None:
            beta1prime_arr = np.array(self.params[4])

        if beta2prime_arr is None:
            beta2prime_arr = np.array(self.params[5])
        params = self.params

        peak1 = np.array([])
        peak2 = np.array([])
        if fixed == 1:
            for beta1prime in beta1prime_arr:
                params[4] = beta1prime
                for beta2prime in beta2prime_arr:
                    params[5] = beta2prime
                    siir = mdl.SIIR(self.SIIR0, params, self.t_sim)
                    siir.runSim()
                    peak1 = np.append(peak1, siir.getDisease1()[1])
                    peak2 = np.append(peak2, siir.getDisease2()[1])
        elif fixed == 2:
            for beta2prime in beta2prime_arr:
                params[5] = beta2prime
                for beta1prime in beta1prime_arr:
                    params[4] = beta1prime
                    siir = mdl.SIIR(self.SIIR0, params, self.t_sim)
                    siir.runSim()
                    peak1 = np.append(peak1, siir.getDisease1()[1])
                    peak2 = np.append(peak2, siir.getDisease2()[1])

        # Plot
        fontsize = 14
        fig, ax = plt.subplots(1,2)
        ax[0].set_xlabel('Beta1\'')
        ax[0].set_ylabel('Beta2\'')
        ax[0].set_title('Disease 1')
        ax[0].set_aspect('equal')
        ax[0].plot(self.params[0] * np.ones(len(beta1prime_arr)), beta2prime_arr, '-.r')
        ax[0].plot(beta1prime_arr, self.params[1] * np.ones(len(beta2prime_arr)), '-r')
        ax[0].legend(('Beta1', 'Beta2'),fontsize=fontsize)
        Beta1prime, Beta2prime = np.meshgrid(beta1prime_arr, beta2prime_arr)
        cf1 = ax[0].contourf(Beta1prime, Beta2prime, peak1, 30)
        ax[0].contour(Beta1prime, Beta2prime, peak1, 30, colors='black')
        fig.colorbar(cf1, ax=ax[0])
        txt = ("Beta1 = " + str(self.params[0]) + "\n Delta1 = " + str(self.params[2]))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        ax[0].text(0.6, 0.4, txt, ha='left', transform=ax[0].transAxes, fontsize=fontsize, bbox=props, wrap=True)
        plt.show()


def testAnalysis():
    # Originals values
    N = 21739040
    SIIR0 = np.zeros(9)
    SIIR0[1] = 119
    SIIR0[2] = 7
    SIIR0[0] = N - np.sum(SIIR0[1:8])

    t_start = 0
    t_end = 15
    n_int = 10000

    t_sim = np.linspace(t_start, t_end, n_int)
    params = [14.3553504, 206.18920775, 13.16275817, 206.04507902,
              14.51848543, 222.29093379, 13.98388944, 154.4459061]

    siirSim_orig = modelAnalysis(SIIR0, params, t_sim)
    siirSim_orig.plotParamsChange()

testAnalysis()