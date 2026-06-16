import numpy as np
import matplotlib.pyplot as plt

fig, axa = plt.subplots()
fig, axc = plt.subplots()
for sid in range(1, 3):
    data = np.load(f"field-aberrations_s{sid}z_prj.pkl", allow_pickle=True)
    z = np.asarray([x[0] for x in data])
    astigmatism = np.asarray([np.hypot(x[1][4], x[1][5]) for x in data])
    astigmatism_5 = np.asarray([x[1][4] for x in data])
    astigmatism_6 = np.asarray([x[1][5] for x in data])
    p_astigmatism = np.polyfit(z, astigmatism, 2)
    g = 8.365 / 25.5
    w222 = p_astigmatism[0] * ((180 * 60 / np.pi) ** 2 * 1e-9 * 2 * np.sqrt(6)) / g**2
    print(f"S{sid}: w222={w222}")

    axa.plot(z, astigmatism_5,  label=f"S{sid}")

    coma = np.asarray([x[1][7] for x in data] )
    axc.plot(z, coma,  label=f"S{sid}")


    # ax.plot(z, np.polyval(p_astigmatism, z), "C3--")
    # ax.plot(z,c222*z**2,'C3--')

axa.set_xlabel("Field-angle [']")
axa.set_ylabel("astigmatism 5 [nm]")
axa.grid()
axa.legend()
axc.set_xlabel("Field-angle [']")
axc.set_ylabel("coma [nm]")
axc.grid()
axc.legend()

