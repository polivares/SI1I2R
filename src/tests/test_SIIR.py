__author__ = "Patricio Andrés Olivares Roncagliolo"
__email__ = "patricio.olivaresr@usm.cl"

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