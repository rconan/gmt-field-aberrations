import os
import numpy as np
import matplotlib.pyplot as plt

home = os.environ["HOME"]
root = home + "/mnt/gmt-point-field-display/"

filename = "gmt-point-field-display_M2S1_z3z6z9.pkl"
print(filename)
probes = np.load(root + filename, allow_pickle=True)

a5 = np.array(
    [
        [x["coefficients"][4] for x in w]
        for w in [x for x in [y["projections"] for y in probes]]
    ]
)
a6 = np.array(
    [
        [x["coefficients"][5] for x in w]
        for w in [x for x in [y["projections"] for y in probes]]
    ]
)

r = np.arange(-100, 101, 20)
fig, ax = plt.subplots()
ax.plot(r, a5)
ax.legend(range(1, 4))
ax.set_prop_cycle(None)
ax.plot(r, a6, "--")
ax.grid()
