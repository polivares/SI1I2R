__author__ = "Patricio Andrés Olivares Roncagliolo"
__email__ = "patricio.olivaresr@usm.cl"

import SIRmodels as mdl
import numpy as np
import matplotlib.pyplot as plt

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
siirSim_orig.runSim()

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
ax[0].legend(('Original disease','Smaller beta2prime', 'Bigger beta2prime'))

ax[1].plot(t_sim, siirSim_orig.getDisease2()[0], '-r')
ax[1].plot(t_sim, siirSim_1.getDisease2()[0], '-g')
ax[1].plot(t_sim, siirSim_2.getDisease2()[0], '-b')
ax[1].legend(('Original disease', 'Smaller beta2prime', 'Bigger beta2prime'))

plt.suptitle("Changes disease 2")

plt.show()