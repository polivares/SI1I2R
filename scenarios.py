__author__ = "Patricio Andrés Olivares Roncagliolo"
__email__ = "patricio.olivaresr@usm.cl"

import SIRmodels as mdl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('TkAgg')

def realScenario():
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

    siirSim_orig = mdl.SIIR(SIIR0, params, t_sim)
    siirSim_orig.runEvaluation()

    # Changes Disease 1
    params = [14.3553504, 206.18920775, 13.16275817, 206.04507902,
              0, 222.29093379, 13.98388944, 154.4459061]

    siirSim_1 = mdl.SIIR(SIIR0, params, t_sim)
    siirSim_1.runEvaluation()

    params = [14.3553504, 206.18920775, 13.16275817, 206.04507902,
              100, 222.29093379, 13.98388944, 154.4459061]

    siirSim_2 = mdl.SIIR(SIIR0, params, t_sim)
    siirSim_2.runEvaluation()

    fig, ax = plt.subplots(1, 2)

    ax[0].plot(t_sim, siirSim_orig.getDisease1()[0], '-r')
    ax[0].plot(t_sim, siirSim_1.getDisease1()[0], '-g')
    ax[0].plot(t_sim, siirSim_2.getDisease1()[0], '-b')
    ax[0].legend(('Original disease', 'Smaller beta1prime', 'Bigger beta1prime'))

    ax[1].plot(t_sim, siirSim_orig.getDisease2()[0], '-r')
    ax[1].plot(t_sim, siirSim_1.getDisease2()[0], '-g')
    ax[1].plot(t_sim, siirSim_2.getDisease2()[0], '-b')
    ax[1].legend(('Original disease', 'Smaller beta1prime', 'Bigger beta1prime'))

    plt.suptitle("Changes disease 1")

    plt.show()

    # Changes Disease 2
    params = [14.3553504, 206.18920775, 13.16275817, 206.04507902,
              14.51848543, 50, 13.98388944, 154.4459061]

    siirSim_1 = mdl.SIIR(SIIR0, params, t_sim)
    siirSim_1.runEvaluation()

    params = [14.3553504, 206.18920775, 13.16275817, 206.04507902,
              14.51848543, 800, 13.98388944, 154.4459061]

    siirSim_2 = mdl.SIIR(SIIR0, params, t_sim)
    siirSim_2.runEvaluation()

    fig, ax = plt.subplots(1, 2)

    ax[0].plot(t_sim, siirSim_orig.getDisease1()[0], '-r')
    ax[0].plot(t_sim, siirSim_1.getDisease1()[0], '-g')
    ax[0].plot(t_sim, siirSim_2.getDisease1()[0], '-b')
    ax[0].legend(('Original disease', 'Smaller beta2prime', 'Bigger beta2prime'))

    ax[1].plot(t_sim, siirSim_orig.getDisease2()[0], '-r')
    ax[1].plot(t_sim, siirSim_1.getDisease2()[0], '-g')
    ax[1].plot(t_sim, siirSim_2.getDisease2()[0], '-b')
    ax[1].legend(('Original disease', 'Smaller beta2prime', 'Bigger beta2prime'))

    plt.suptitle("Changes disease 2")

    plt.show()

def getScenario():
    # Originals values
    N = 1000000
    SIIR0 = np.zeros(9)
    SIIR0[1] = 100
    SIIR0[2] = 7
    SIIR0[0] = N - np.sum(SIIR0[1:8])

    t_start = 0
    t_end = 50
    n_int = 10000

    t_sim = np.linspace(t_start, t_end, n_int)
    params = [10, 10, 5, 9.9, 10, 10, 3, 8]

    siirSim_orig = mdl.SIIR(SIIR0, params, t_sim)
    siirSim_orig.runEvaluation()

    # Changes Disease 1
    ###################

    # Smaller
    params = [10, 10, 5, 9.9, 2, 10, 3, 8]

    siirSim_1 = mdl.SIIR(SIIR0, params, t_sim)
    siirSim_1.runEvaluation()

    # Bigger
    params = [10, 10, 5, 9.9, 20, 10, 3, 8]

    siirSim_2 = mdl.SIIR(SIIR0, params, t_sim)
    siirSim_2.runEvaluation()

    fig, ax = plt.subplots(1, 2)

    ax[0].plot(t_sim, siirSim_orig.getDisease1()[0], '-r')
    ax[0].plot(t_sim, siirSim_1.getDisease1()[0], '-g')
    ax[0].plot(t_sim, siirSim_2.getDisease1()[0], '-b')
    ax[0].legend(('Original disease', 'Smaller beta1prime', 'Bigger beta1prime'))

    ax[1].plot(t_sim, siirSim_orig.getDisease2()[0], '-r')
    ax[1].plot(t_sim, siirSim_1.getDisease2()[0], '-g')
    ax[1].plot(t_sim, siirSim_2.getDisease2()[0], '-b')
    ax[1].legend(('Original disease', 'Smaller beta1prime', 'Bigger beta1prime'))

    plt.suptitle("Changes disease 1")

    plt.show()

    # Changes Disease 2
    ###################

    # Smaller
    params = [10, 10, 5, 9.9, 10, 2, 3, 8]

    siirSim_1 = mdl.SIIR(SIIR0, params, t_sim)
    siirSim_1.runEvaluation()

    # Bigger
    params = [10, 10, 5, 9.9, 10, 20, 3, 8]

    siirSim_2 = mdl.SIIR(SIIR0, params, t_sim)
    siirSim_2.runEvaluation()

    fig, ax = plt.subplots(1, 2)

    ax[0].plot(t_sim, siirSim_orig.getDisease1()[0], '-r')
    ax[0].plot(t_sim, siirSim_1.getDisease1()[0], '-g')
    ax[0].plot(t_sim, siirSim_2.getDisease1()[0], '-b')
    ax[0].legend(('Original disease', 'Smaller beta2prime', 'Bigger beta2prime'))

    ax[1].plot(t_sim, siirSim_orig.getDisease2()[0], '-r')
    ax[1].plot(t_sim, siirSim_1.getDisease2()[0], '-g')
    ax[1].plot(t_sim, siirSim_2.getDisease2()[0], '-b')
    ax[1].legend(('Original disease', 'Smaller beta2prime', 'Bigger beta2prime'))

    plt.suptitle("Changes disease 2")

    #plt.show()

def eval1():
    # Originals values
    N = 1000000
    SIIR0 = np.zeros(9)
    SIIR0[1] = 10
    SIIR0[2] = 10
    SIIR0[0] = N - np.sum(SIIR0[1:8])

    t_start = 0
    t_end = 1.3
    n_int = 10000

    t_sim = np.linspace(t_start, t_end, n_int)
    #params = [10, 10, 8, 8.7, 10, 10, 8, 8.7]
    params = [100, 100, 80, 87, 100, 100, 80, 87]

    siirSim_orig = mdl.SIIR(SIIR0, params, t_sim)
    siirSim_orig.runEvaluation()

    # Smaller
    #params = [10, 10, 8, 8.7, 40, 40, 8, 8.7]
    params = [100, 100, 80, 87, 400, 400, 80, 87]

    siirSim_1 = mdl.SIIR(SIIR0, params, t_sim)
    siirSim_1.runEvaluation()

    fig, ax = plt.subplots(1, 1)
    #ax.set_yticklabels([])
    #ax.set_xticklabels([])
    #ax.tick_params(axis=u'both', which=u'both',length=0)
    ax.plot(t_sim, siirSim_orig.getDisease1()[0]/N, '-r')
    ax.plot(t_sim, siirSim_1.getDisease1()[0]/N, '-.r')
    ax.plot(t_sim, siirSim_orig.getDisease2()[0]/N, '-b')
    ax.plot(t_sim, siirSim_1.getDisease2()[0]/N, '-.b')

    #ax.annotate("", xy=(0.5, 0.5), xytext=(0, 0), arrowprops=dict(arrowstyle="->"))
    #ax.arrow(5.21, 19311, 4.66-5.21, 23388-19311,head_width=0.14, head_length=600, length_includes_head = True)

    ax.legend(('Disease1', 'Disease 1 with interaction', 'Disease2', 'Disease2 with interaction'))
    plt.title("Cooperative Interaction")
    #plt.show()
    plt.savefig("/home/polivares/cooperative.png", dpi=900)

def eval2():
    # Originals values
    N = 1000000
    SIIR0 = np.zeros(9)
    SIIR0[1] = 10
    SIIR0[2] = 10
    SIIR0[0] = N - np.sum(SIIR0[1:8])

    t_start = 0
    t_end = 1.3
    n_int = 10000

    t_sim = np.linspace(t_start, t_end, n_int)
    params = [100, 100, 80, 87, 100, 100, 80, 87]

    siirSim_orig = mdl.SIIR(SIIR0, params, t_sim)
    siirSim_orig.runEvaluation()

    # Smaller
    params = [100, 100, 80, 87, 0, 0, 80, 87]

    siirSim_1 = mdl.SIIR(SIIR0, params, t_sim)
    siirSim_1.runEvaluation()

    fig, ax = plt.subplots(1, 1)
    #ax.set_yticklabels([])
    #ax.set_xticklabels([])
    #ax.tick_params(axis=u'both', which=u'both', length=0)
    ax.plot(t_sim, siirSim_orig.getDisease1()[0]/N, '-r')
    ax.plot(t_sim, siirSim_1.getDisease1()[0]/N, '-.r')
    ax.plot(t_sim, siirSim_orig.getDisease2()[0]/N, '-b')
    ax.plot(t_sim, siirSim_1.getDisease2()[0]/N, '-.b')
    #ax.annotate("", xy=(0.5, 0.5), xytext=(0, 0), arrowprops=dict(arrowstyle="->"))
    #ax.arrow(5.21, 19311, 4.66-5.21, 23388-19311,head_width=0.14, head_length=600, length_includes_head = True)

    ax.legend(('Disease1', 'Disease 1 with interaction', 'Disease2', 'Disease2 with interaction'))
    plt.title("Competitive interaction")
    #plt.show()
    plt.savefig("/home/polivares/coompetitive.png", dpi=900)

eval1()
eval2()