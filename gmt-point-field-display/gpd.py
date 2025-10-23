import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

home = os.environ["HOME"]
root = home + "/mnt/gmt-point-field-display/M2S1_tx/"

for i, tx in enumerate(range(-100, 101, 20)):

    filename = "gmt-point_field-display.pkl"
    print(filename)
    probes = np.load(root + filename, allow_pickle=True)

    z = np.asarray([x["pointing"]["zenith"]["Arcminute"] for x in probes["fields"]])
    a = np.asarray([x["pointing"]["azimuth"]["Radian"] for x in probes["fields"]])
    x = z * np.cos(a)
    y = z * np.sin(a)
    triang = tri.Triangulation(x, y)

    z5 = np.asarray([x["coefficients"][4] for x in probes["projections"]])
    z6 = np.asarray([x["coefficients"][5] for x in probes["projections"]])
    z56 = np.hypot(z5, z6)

    fig, ax = plt.subplots()
    h = ax.tripcolor(triang, z56, shading="gouraud", cmap="Spectral")
    ax.set_aspect("equal")
    ax.grid()
    ax.set_title(f"M2 S1 - Tx={tx}micron")
    fig.savefig(root + f"gffd{i:02d}.png", bbox_inches="tight")
    plt.close(fig)
