import matplotlib.pyplot as plt
import pandas as pd

fig, ax = plt.subplots(1, 2, figsize=(10, 5))

T = 300
isotropic = pd.read_csv("isotropic_T_{}.csv".format(T))
simohughes = pd.read_csv("simohughes_T_{}.csv".format(T))
ax[0].plot(isotropic["strain"], isotropic["stress"], label="old tensor mechanics")
ax[0].plot(simohughes["strain"], simohughes["stress"], label="new tensor mechanics")
ax[0].set_xlabel("strain")
ax[0].set_ylabel("stress")
ax[0].legend()
ax[0].set_title("T = {}".format(T))

T = 2000
isotropic = pd.read_csv("isotropic_T_{}.csv".format(T))
simohughes = pd.read_csv("simohughes_T_{}.csv".format(T))
ax[1].plot(isotropic["strain"], isotropic["stress"], label="old tensor mechanics")
ax[1].plot(simohughes["strain"], simohughes["stress"], label="new tensor mechanics")
ax[1].set_xlabel("strain")
ax[1].set_ylabel("stress")
ax[1].legend()
ax[1].set_title("T = {}".format(T))

plt.show()
