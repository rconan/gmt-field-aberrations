import numpy as np
import matplotlib.pyplot as plt
from matplotlib import tri

file = "zernike-modes_5_gs0_s2_z10.0a0.0.pkl"
data = np.load(file, allow_pickle=True)

xy = data["xy"]
x = np.array([c[0] for c in xy])
y = np.array([c[1] for c in xy])
modes = np.array(data["modes"]).reshape(-1, len(xy))
min = modes.min()
max = modes.max()

mesh = tri.Triangulation(x, y)
fig, ax = plt.subplots()
k = 0
for i in range(5):
    for j in range(i):
        x0 = j
        y0 = i
        h = ax.tripcolor(
            x0 + mesh.x,
            y0 + mesh.y,
            modes[k, :],
            vmin=min,
            vmax=max,
            cmap="nipy_spectral",
        )
        k += 1
ax.set_title(file)
ax.set_aspect("equal")
fig.colorbar(h, ax=ax)
