import numpy as np
import matplotlib.pyplot as plt

opdset = np.load("opds.pkl", allow_pickle=True)

delta_max = 1e9*np.max([np.array([opd["delta"] for opd in opds]).max() for opds in opdset])
delta_min = 1e9*np.min([np.array([opd["delta"] for opd in opds]).min() for opds in opdset])

fig, ax = plt.subplots()
for opds in opdset:
    x = np.array([opd["xyz"][0] for opd in opds])
    y = np.array([opd["xyz"][1] for opd in opds])
    delta = np.array([opd["delta"] for opd in opds])
    h = ax.tripcolor(x, y, delta * 1e9, vmin=delta_min, vmax=delta_max)
    ax.set_aspect("equal")

ax.grid()
ax.set_aspect("equal")
fig.colorbar(h, ax=ax)
