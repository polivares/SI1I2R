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


    def plotPeak2D(self, beta1prime_arr=None, beta2prime_arr=None, savefig=False):
        if beta1prime_arr is None:
            beta1prime_arr = np.array([self.params[4]])

        if beta2prime_arr is None:
            beta2prime_arr = np.array([self.params[5]])
        params = self.params.copy()

        Beta1prime, Beta2prime = np.meshgrid(beta2prime_arr, beta1prime_arr)
        peak1 = np.zeros((len(beta1prime_arr), len(beta2prime_arr)))
        peak2 = np.zeros((len(beta1prime_arr), len(beta2prime_arr)))
        for j in np.arange(len(beta2prime_arr)):
            beta2prime = beta2prime_arr[j]
            params[5] = beta2prime
            for i in np.arange(len(beta1prime_arr)):
                beta1prime = beta1prime_arr[i]
                print("(beta1prime,beta2prime)", beta1prime, beta2prime)
                params[4] = beta1prime
                siir = mdl.SIIR(self.SIIR0, params, self.t_sim)
                siir.runEvaluation()
                peak1[i, j] = self.t_sim[siir.getDisease1()[1]][0] # getPeak1
                peak2[i, j] = self.t_sim[siir.getDisease2()[1]][0] # getPeak2

        # Plot
        fontsize = 8
        fig, ax = plt.subplots(1, 2)
        ax[0].set_xlabel('Beta1\'')
        ax[0].set_ylabel('Beta2\'')
        ax[0].set_title('Disease 1')
        ax[0].set_aspect('equal')
        ax[0].plot(self.params[0] * np.ones(len(beta2prime_arr)), beta2prime_arr, '-.r')
        ax[0].plot(beta1prime_arr, self.params[1] * np.ones(len(beta1prime_arr)), '-r')
        ax[0].legend(('Beta1', 'Beta2'), fontsize=fontsize)
        cf1 = ax[0].contourf(Beta2prime, Beta1prime, peak1, 30, cmap='hot')
        #ax[0].contour(Beta1prime, Beta2prime, peak1, 10, colors='black')
        fig.colorbar(cf1, ax=ax[0])
        txt = ("Beta1 = " + str(self.params[0]) + "\n Delta1 = " + str(self.params[2])
               + "\n Beta1\' = " + str(self.params[4]) + "\n Delta1\' = " + str(self.params[6]))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        ax[0].text(0.6, 0.4, txt, ha='left', transform=ax[0].transAxes, fontsize=fontsize, bbox=props, wrap=True)

        ax[1].set_xlabel('Beta1\'')
        ax[1].set_ylabel('Beta2\'')
        ax[1].set_title('Disease 2')
        ax[1].set_aspect('equal')
        ax[1].plot(self.params[0] * np.ones(len(beta2prime_arr)), beta2prime_arr, '-.r')
        ax[1].plot(beta1prime_arr, self.params[1] * np.ones(len(beta1prime_arr)), '-r')
        ax[1].legend(('Beta1', 'Beta2'), fontsize=fontsize)
        cf2 = ax[1].contourf(Beta2prime, Beta1prime, peak2, 30, cmap='winter')
        #ax[1].contour(Beta1prime, Beta2prime, peak2, 10, colors='black')
        fig.colorbar(cf2, ax=ax[1])
        txt = ("Beta2 = " + str(self.params[1]) + "\n Delta2 = " + str(self.params[3])
               + "\n Beta2\' = " + str(self.params[5]) + "\n Delta2\' = " + str(self.params[7]))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        ax[1].text(0.6, 0.4, txt, ha='left', transform=ax[1].transAxes, fontsize=fontsize, bbox=props, wrap=True)

        if savefig:
            fig.savefig("images/peak2D.svg", format='svg')
        else:
            plt.show()

    def plotInfected2D(self, beta1prime_arr=None, beta2prime_arr=None, savefig=False):
        if beta1prime_arr is None:
            beta1prime_arr = np.array([self.params[4]])

        if beta2prime_arr is None:
            beta2prime_arr = np.array([self.params[5]])
        params = self.params.copy()

        Beta1prime, Beta2prime = np.meshgrid(beta2prime_arr, beta1prime_arr)
        inf1 = np.zeros((len(beta2prime_arr), len(beta1prime_arr)))
        inf2 = np.zeros((len(beta2prime_arr), len(beta1prime_arr)))
        for j in np.arange(len(beta2prime_arr)):
            beta2prime = beta2prime_arr[j]
            params[5] = beta2prime
            for i in np.arange(len(beta1prime_arr)):
                beta1prime = beta1prime_arr[i]
                print("(beta1prime,beta2prime)", beta1prime, beta2prime)
                params[4] = beta1prime
                siir = mdl.SIIR(self.SIIR0, params, self.t_sim)
                siir.runEvaluation()
                inf1[i, j] = siir.getNInfected1() # n infected 1
                inf2[i, j] = siir.getNInfected2() # n infected 2

        # Plot
        fontsize = 8
        fig, ax = plt.subplots(1, 2)
        ax[0].set_xlabel('Beta1\'')
        ax[0].set_ylabel('Beta2\'')
        ax[0].set_title('Disease 1')
        ax[0].set_aspect('equal')
        ax[0].plot(self.params[0] * np.ones(len(beta2prime_arr)), beta2prime_arr, '-.r')
        ax[0].plot(beta1prime_arr, self.params[1] * np.ones(len(beta1prime_arr)), '-r')
        ax[0].legend(('Beta1', 'Beta2'), fontsize=fontsize)
        cf1 = ax[0].contourf(Beta2prime, Beta1prime, inf1, 30, cmap='hot')
        #ax[0].contour(Beta1prime, Beta2prime, peak1, 10, colors='black')
        fig.colorbar(cf1, ax=ax[0])
        txt = ("Beta1 = " + str(self.params[0]) + "\n Delta1 = " + str(self.params[2])
               + "\n Beta1\' = " + str(self.params[4]) + "\n Delta1\' = " + str(self.params[6]))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        ax[0].text(0.6, 0.4, txt, ha='left', transform=ax[0].transAxes, fontsize=fontsize, bbox=props, wrap=True)

        ax[1].set_xlabel('Beta1\'')
        ax[1].set_ylabel('Beta2\'')
        ax[1].set_title('Disease 2')
        ax[1].set_aspect('equal')
        ax[1].plot(self.params[0] * np.ones(len(beta2prime_arr)), beta2prime_arr, '-.r')
        ax[1].plot(beta1prime_arr, self.params[1] * np.ones(len(beta1prime_arr)), '-r')
        ax[1].legend(('Beta1', 'Beta2'), fontsize=fontsize)
        cf2 = ax[1].contourf(Beta2prime, Beta1prime, inf2, 30, cmap='winter')
        #ax[1].contour(Beta1prime, Beta2prime, peak2, 10, colors='black')
        fig.colorbar(cf2, ax=ax[1])
        txt = ("Beta2 = " + str(self.params[1]) + "\n Delta2 = " + str(self.params[3])
               + "\n Beta2\' = " + str(self.params[5]) + "\n Delta2\' = " + str(self.params[7]))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        ax[1].text(0.6, 0.4, txt, ha='left', transform=ax[1].transAxes, fontsize=fontsize, bbox=props, wrap=True)

        if savefig:
            fig.savefig("images/infected2D.svg", format='svg')
        else:
            plt.show()

    def plotPeakNInfected(self, beta1prime_arr=None, beta2prime_arr=None, savefig=False):
        if beta1prime_arr is None:
            beta1prime_arr = np.array([self.params[4]])

        if beta2prime_arr is None:
            beta2prime_arr = np.array([self.params[5]])
        params = self.params.copy()

        Beta1prime, Beta2prime = np.meshgrid(beta2prime_arr, beta1prime_arr)
        peak1 = np.zeros((len(beta1prime_arr), len(beta2prime_arr)))
        peak2 = np.zeros((len(beta1prime_arr), len(beta2prime_arr)))
        inf1 = np.zeros((len(beta2prime_arr), len(beta1prime_arr)))
        inf2 = np.zeros((len(beta2prime_arr), len(beta1prime_arr)))
        for j in np.arange(len(beta2prime_arr)):
            beta2prime = beta2prime_arr[j]
            params[5] = beta2prime
            for i in np.arange(len(beta1prime_arr)):
                beta1prime = beta1prime_arr[i]
                print("(beta1prime,beta2prime)", beta1prime, beta2prime)
                params[4] = beta1prime
                siir = mdl.SIIR(self.SIIR0, params, self.t_sim)
                siir.runEvaluation()
                peak1[i, j] = self.t_sim[siir.getDisease1()[1]][0] # getPeak1
                peak2[i, j] = self.t_sim[siir.getDisease2()[1]][0] # getPeak2
                inf1[i, j] = siir.getNInfected1()  # n infected 1
                inf2[i, j] = siir.getNInfected2()  # n infected 2
        peak1 = peak1.flatten()
        peak2 = peak2.flatten()
        inf1 = inf1.flatten()
        inf2 = inf2.flatten()

        plt.plot(peak1, inf1, '.r', label='Disease 1')
        plt.xlabel("Peak position")
        plt.ylabel("N° people infected")
        plt.show()

        plt.plot(peak2, inf2, '.b', label='Disease 2')
        plt.xlabel("Peak position")
        plt.ylabel("N° people infected")
        plt.show()



    def plotChangeBoth(self, beta1prime_arr=None, beta2prime_arr=None):
        if beta1prime_arr is None:
            beta1prime_arr = np.array([self.params[4]])

        if beta2prime_arr is None:
            beta2prime_arr = np.array([self.params[5]])


        fig,ax=plt.subplots(1,2)
        siir = mdl.SIIR(self.SIIR0, self.params, self.t_sim)
        siir.runEvaluation()
        ax[0].plot(self.t_sim, siir.getDisease1()[0], '-r', label='Original')
        ax[1].plot(self.t_sim, siir.getDisease2()[0], '-b', label='Original')

        params = self.params.copy()
        for beta1prime in beta1prime_arr[::-1]:
            for beta2prime in beta2prime_arr[::-1]:
                params[4] = beta1prime
                params[5] = beta2prime
                siir = mdl.SIIR(self.SIIR0, params, self.t_sim)
                siir.runEvaluation()
                ax[0].plot(self.t_sim, siir.getDisease1()[0], label='B1\'=' + str(beta1prime) + "\nB2\'=" + str(beta2prime))
                ax[1].plot(self.t_sim, siir.getDisease2()[0], label='B1\'=' + str(beta1prime) + "\nB2\'=" + str(beta2prime))
        ax[0].legend()
        ax[1].legend()
        txt = ("Beta1 = " + str(self.params[0]) + ", Delta1 = " + str(self.params[2]) +
               ", Beta2 = " + str(self.params[1]) + ", Delta2 = " + str(self.params[3]) +
               "\nBeta1' = " + str(self.params[4]) + ", Delta1' = " + str(self.params[6]) +
               ", Beta2' = " + str(self.params[5]) + ", Delta2' = " + str(self.params[7]))
        plt.figtext(0.65, 0.9, txt, ha="right", va="bottom", bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        plt.show()

    def plotChange(self, betaprime_arr=None, fixed=2, peak=False):
        if betaprime_arr is None:
            if fixed == 1:
                betaprime_arr = np.array([self.params[5]])
            elif fixed == 2:
                betaprime_arr = np.array([self.params[4]])

        fig, ax = plt.subplots(2, 1)
        siir = mdl.SIIR(self.SIIR0, self.params, self.t_sim)
        siir.runEvaluation()
        ax[0].plot(self.t_sim, siir.getDisease1()[0], '-r', label='Original')
        ax[1].plot(self.t_sim, siir.getDisease2()[0], '-r', label='Original')

        if peak:
            ax[0].plot(self.t_sim[siir.getDisease1()[1]], siir.getDisease1()[0][siir.getDisease1()[1]], '*k')
            ax[1].plot(self.t_sim[siir.getDisease2()[1]], siir.getDisease2()[0][siir.getDisease2()[1]], '*k')

        params = self.params.copy()
        if fixed == 2:
            for beta1prime in betaprime_arr[::-1]:
                params[4] = beta1prime
                siir = mdl.SIIR(self.SIIR0, params, self.t_sim)
                siir.runEvaluation()
                ax[0].plot(self.t_sim, siir.getDisease1()[0], label='Beta1\'=' + str(beta1prime))
                ax[1].plot(self.t_sim, siir.getDisease2()[0], label='Beta1\'=' + str(beta1prime))
                if peak:
                    ax[0].plot(self.t_sim[siir.getDisease1()[1]], siir.getDisease1()[0][siir.getDisease1()[1]], '*k')
                    ax[1].plot(self.t_sim[siir.getDisease2()[1]], siir.getDisease2()[0][siir.getDisease2()[1]], '*k')

        elif fixed == 1:
            for beta2prime in betaprime_arr[::-1]:
                params[5] = beta2prime
                siir = mdl.SIIR(self.SIIR0, params, self.t_sim)
                siir.runEvaluation()
                ax[1].plot(self.t_sim, siir.getDisease2()[0], label='Beta2\'=' + str(beta2prime))
                ax[0].plot(self.t_sim, siir.getDisease1()[0], label='Beta2\'=' + str(beta2prime))
                if peak:
                    ax[0].plot(self.t_sim[siir.getDisease1()[1]], siir.getDisease1()[0][siir.getDisease1()[1]], '*k')
                    ax[1].plot(self.t_sim[siir.getDisease2()[1]], siir.getDisease2()[0][siir.getDisease2()[1]], '*k')


        ax[0].legend()
        ax[1].legend()
        txt = ("Beta1 = " + str(self.params[0]) + ", Delta1 = " + str(self.params[2]) +
               ", Beta2 = " + str(self.params[1]) + ", Delta2 = " + str(self.params[3]) +
               "\nBeta1' = " + str(self.params[4]) + ", Delta1' = " + str(self.params[6]) +
               ", Beta2' = " + str(self.params[5]) + ", Delta2' = " + str(self.params[7]))
        plt.figtext(0.65, 0.9, txt, ha="right", va="bottom",
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
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
    # siirSim_orig.plotChange2D(beta1prime_arr=np.arange(0, 40, 0.1), beta2prime_arr=np.arange(180, 250, 1))
    siirSim_orig.plotChange(betaprime_arr=np.arange(200, 240, 5), fixed=1)
    siirSim_orig.plotChange(betaprime_arr=np.arange(5, 15, 1), fixed=2)

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

    params = [10, 10, 9, 5, 10, 10, 15, 5] # Fixed delta primes
    #params = [10, 10, 5, 9, 10, 10, 5, 5]  # Fixed delta primes


    siirSim_orig = modelAnalysis(SIIR0, params, t_sim)

    #siirSim_orig.plotChange(betaprime_arr=np.array([5, 15, 40]), fixed=1, peak=False)
    #siirSim_orig.plotChange(betaprime_arr=np.array([5, 15, 40]), fixed=2, peak=False)
    #siirSim_orig.plotPeak2D(beta1prime_arr=np.arange(0, 40, 2), beta2prime_arr=np.arange(0, 40, 2), savefig=False)
    #siirSim_orig.plotInfected2D(beta1prime_arr=np.arange(0, 40, 2), beta2prime_arr=np.arange(0, 40, 2), savefig=False)
    siirSim_orig.plotPeakNInfected(beta1prime_arr=np.arange(0, 40, 2), beta2prime_arr=np.arange(0, 40, 2), savefig=False)

getAnalysis()