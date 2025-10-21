#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
# %matplotlib inline

probes = np.load("gmt-full-field.pkl", allow_pickle=True)

z = np.asarray([x["pointing"]["zenith"]["Arcminute"] for x in probes["fields"]])
a = np.asarray([x["pointing"]["azimuth"]["Radian"] for x in probes["fields"]])
x = z * np.cos(a)
y = z * np.sin(a)
triang = tri.Triangulation(x, y)

fig, ax = plt.subplots()
ax.triplot(triang)
ax.set_aspect("equal")
ax.grid()

z5 = np.asarray([x["coefficients"][4] for x in probes["projections"]])
z6 = np.asarray([x["coefficients"][5] for x in probes["projections"]])
z7 = np.asarray([x["coefficients"][6] for x in probes["projections"]])
z8 = np.asarray([x["coefficients"][7] for x in probes["projections"]])

fig, axs = plt.subplots(ncols=2, nrows=2, sharex=True, sharey=True, figsize=(12, 10))
ax = axs[0, 0]
h = ax.tripcolor(triang, z5, shading="gouraud", cmap="Spectral")
ax.set_aspect("equal")
ax.grid()
fig.colorbar(h, ax=ax)
ax.set_title("Z5")
ax = axs[0, 1]
h = ax.tripcolor(triang, z6, shading="gouraud", cmap="Spectral")
ax.set_aspect("equal")
ax.grid()
fig.colorbar(h, ax=ax)
ax.set_title("Z6")
ax = axs[1, 0]
h = ax.tripcolor(triang, z7, shading="gouraud", cmap="Spectral")
ax.set_aspect("equal")
ax.grid()
fig.colorbar(h, ax=ax)
ax.set_title("Z7")
ax = axs[1, 1]
h = ax.tripcolor(triang, z8, shading="gouraud", cmap="Spectral")
ax.set_aspect("equal")
ax.grid()
fig.colorbar(h, ax=ax)
ax.set_title("Z8")

probes["fields"][0]["pupil_mode"]


# In[ ]:
