import numpy as np
import matplotlib.pyplot as plt

field_zerns = np.load("zerns.pkl", allow_pickle=True)

zen = np.asarray([fz["za"][0] for fz in field_zerns])
azi = np.asarray([fz["za"][1] for fz in field_zerns])
(x, y) = (zen * np.cos(azi), zen * np.sin(azi))
modes = [fz["modes"] for fz in field_zerns]
n_mode = len(modes[0])
coefs = [[m["coef"] for m in mode] for mode in modes]


fig, axs = plt.subplots(ncols=n_mode,figsize=(12,5))
for i, ax in enumerate(axs):
    c = [c[i] for c in coefs]
    h = ax.tripcolor(x, y, c, shading="gouraud")
    ax.grid()
    ax.set_aspect("equal")
    fig.colorbar(h, ax=ax)
