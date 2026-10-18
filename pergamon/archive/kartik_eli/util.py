import numpy as np
import matplotlib.pyplot as plt
import math
#from astropy import constants

'''
Equations taken from:
Rahvar, S., Mehrabi, A., & Dominik, M. 2011, MNRAS, 410, 912
and
Poleski, R., & Yee, J. C. 2019, Astronomy and Computing, 26, 35
'''
# t_star = 0.051
# t_start: t_0 - 5.5 * t_star,
# t_stop: t_0 + 5.5 * t_star,
# subtract_2460000: True,

# speed of light, c, (au/day)
# orbital radius, a, (au)
# mass of source star, m_star, mass of lens, M (M_sun)
# inclination angle with respect to observer-lens line of sight, φ, in radians


def mag_max(M, m_star, a):

    rho_star = 22.7 * ((M) ** -0.5) * ((m_star) ** 0.8) * ((a) ** -0.5)

    mu_max = np.sqrt(1 + 4/(rho_star ** 2))

    return mu_max

M_vals = np.linspace(0.1, 1000.0, 100000)
mstar_vals = np.logspace(-1.0, 3.0, 100000)
orbradius_vals = np.linspace(0.1, 100.0, 100000)

mumax_vals = np.zeros(100000)

def chabrier(m_star):
    if m_star < 1:
        distribution = 0.093 * math.exp( (- ((math.log(m_star) - math.log(0.2)) ** 2)) / (2 * (0.55 ** 2) ))
    elif m_star >= 1:
        distribution = 0.041 * (m_star ** -1.35)

    return distribution

def poleski_main(t, t_0, u_0, t_E):
    # eq5
    tao = (t - t_0)/t_E

    # eq4
    u = np.sqrt((u_0 ** 2) + (tao ** 2))

    # eq5
    Aofu = ((u ** 2) + 2) / (u * np.sqrt((u ** 2) + 4))

    return Aofu

t_01 = 0.0
u_01 = 3.325 # 9.002393432251264
t_E1 = 0.008018544581274847

t_vals1 = np.linspace(-2.5, 2.5, 100000)

amp_vals1 = np.zeros(100000)
amp_vals2 = np.zeros(100000)

def rahvar_main(t, t_0, φ, M, m_star, a):

    c = 173.145 # (au/day)
    G = 0.000295913010 # au * 1/solar mass * (au/day)^2
    R_S = (2 * M * G)/(c ** 2) # au

    # equation 17
    t_E = ((2 * a)/c) * (np.sqrt(M/(m_star + M))) * 1/24
    # equation 18
    tao = (t - t_0)/t_E
    u_0 = np.sqrt(φ * (a / ((2 * R_S))))
    u_t = np.sqrt((u_0 ** 2) + (tao ** 2))

    # equaiton 14
    Amplitude = (2/np.pi) * ((1 + 1/(u_t ** 2)) * np.arcsin(1/(np.sqrt(1 + 1/(u_t ** 2)))) + 1/u_t)

    return Amplitude

t_02 = 0.0
φ = 1.59989 * (10.0 ** (-6)) # radians
M = 8.5
m_star = 0.35
a = 17



additive = 0

for element in t_vals1:
    amp_vals1[additive] = (rahvar_main(element, t_02, φ, M, m_star, a))
    additive += 1

plt.figure(1)
plt.plot(t_vals1, amp_vals1, label="Rahvar")
plt.xlabel("Time [Hours]")
plt.ylabel("Magnfication")


additive = 0

for element in t_vals1:
    amp_vals2[additive] = (poleski_main(element, t_01, u_01, t_E1))
    additive += 1

plt.plot(t_vals1, amp_vals2, label="Poleski")
plt.xlabel("Time [Hours]")
plt.ylabel("Magnification")
plt.yscale("linear")
plt.legend()

plt.figure(2)
residual = (amp_vals1 - amp_vals2) + 1
#print(residual)
#print(amp_vals2)

#print(amp_vals1)

plt.plot(t_vals1, residual, color = "green", label = "Residual")
plt.xlabel("Time [Hours]")
plt.ylabel("Magnification")
plt.legend()
# second set

additive = 0

for element in M_vals:
    mumax_vals[additive] = (mag_max(element, 1, 1))
    additive += 1

plt.figure(3)
plt.plot(M_vals, mumax_vals)
plt.xlabel("Mass of Lens [$M_{\\odot}$]")
plt.ylabel("Maximum Magnification")
plt.title("$M_{*}$ = 1 $M_{\\odot}$, a = 1 AU")
plt.legend()

mumax_vals = np.zeros(100000)

additive = 0

for element in mstar_vals:
    mumax_vals[additive] = (mag_max(1, element, 1))
    additive += 1

plt.figure(4)
plt.plot(mstar_vals, mumax_vals)
plt.xlabel("Mass of Star [$M_{\\odot}$]")
plt.xscale("log")
plt.ylabel("Maximum Magnification")
plt.title("$M$ = 1 $M_{\\odot}$, a = 1 AU")
plt.legend()

mumax_vals = np.zeros(100000)

additive = 0

for element in orbradius_vals:
    mumax_vals[additive] = (mag_max(1, 1, element))
    additive += 1

plt.figure(5)
plt.plot(orbradius_vals, mumax_vals)
plt.xlabel("Orbital Radius [AU]")
#plt.xscale("log")
plt.ylabel("Maximum Magnification")
plt.title("$M$ = 1 $M_{\\odot}$, $M_{*}$ = 1 $M_{\\odot}$")
plt.legend()


orbradius_vals = np.linspace(0.0001, 0.55, 100000)

additive = 0

for element in orbradius_vals:
    mumax_vals[additive] = (mag_max(1, 1, element))
    additive += 1

#rahvar equation 22
# Mstar = 1, M = 1
orbperiod_vals = ((2 * np.pi) * np.sqrt((orbradius_vals ** 3) / (2 * 0.000295913010)))

# orbperiod_vals = np.linspace(0.1, 100.0, 100000)
plt.figure(6)
plt.plot(orbperiod_vals, mumax_vals)
plt.xlabel("Orbital Period [Days]")
#plt.xscale("log")
plt.ylabel("Maximum Magnification")
plt.title("$M$ = 1 $M_{\\odot}$, $M_{*}$ = 1 $M_{\\odot}$")

# chabrier initial mass function, rahvar eq 20

mstar_vals = np.logspace(-2.0, 3.0, 100000)

chabrier_distribution = np.zeros(100000)

additive = 0

for element in mstar_vals:
    chabrier_distribution[additive] = (chabrier(element))
    additive += 1

plt.figure(7)
plt.plot(mstar_vals, chabrier_distribution)
plt.xlabel("Mass of Star [$M_{\\odot}$]")
plt.ylabel("Mass Function")
plt.xscale("log")
plt.yscale("log")
plt.title("Chabrier Distribution")

plt.show()
