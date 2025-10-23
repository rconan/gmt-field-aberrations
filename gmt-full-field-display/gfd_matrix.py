import pickle
import os

home = os.environ["HOME"]
root = home + "/mnt/gmt-full-field-display/matrix/"

if os.path.exists(root + "matrix.pkl"):
    with open(root + "matrix.pkl", "rb") as f:
        data = pickle.load(f)

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri

    z = data["zenith"]
    a = data["azimuth"]
    x = z * np.cos(a)
    y = z * np.sin(a)
    triang = tri.Triangulation(x, y)

    fig, axs = plt.subplots(nrows=4, ncols=7, figsize=(19, 11), sharex=True, sharey=True)
    j = 0
    for r in ["T", "R"]:
        for a in ["x", "y"]:
            for i in range(7):
                z5 = data[f"{r}{a}"][i][0]
                z6 = data[f"{r}{a}"][i][1]
                z56 = np.hypot(z5, z6)
                ax = axs[j, i]
                h = ax.tripcolor(triang, z56, shading="gouraud", cmap="Spectral")
                ax.set_aspect("equal")
                ax.grid()
                if j == 0:
                    ax.set_title(f"S{i + 1}")
                if i == 0:
                    ax.set_ylabel(f"{r}{a}")
            j += 1
    fig.savefig(root + "matrix.png", bbox_inches="tight")
else:
    data = {}
    for r in ["T", "R"]:
        for a in ["x", "y"]:
            data[f"{r}{a}"] = []
            for i in range(1, 8):
                if r == "T":
                    filename = f"gmt-full-field_M2S{i}_T{a}R.pkl"
                else:
                    filename = f"gmt-full-field_M2S{i}_TR{a}.pkl"
                print(filename)
                probes = np.load(root + filename, allow_pickle=True)

                if "zenith" not in data:
                    data["zenith"] = np.asarray(
                        [x["pointing"]["zenith"]["Arcminute"] for x in probes["fields"]]
                    )
                    data["azimuth"] = np.asarray(
                        [x["pointing"]["azimuth"]["Radian"] for x in probes["fields"]]
                    )

                z5 = np.asarray([x["coefficients"][4] for x in probes["projections"]])
                z6 = np.asarray([x["coefficients"][5] for x in probes["projections"]])

                data[f"{r}{a}"] += [(z5, z6)]

    with open(root + "/matrix.pkl", "wb") as f:
        pickle.dump(data, f)
