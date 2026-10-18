import numpy as np
import matplotlib.pyplot as plt

'''
Equations taken from:
Rahvar, S., Mehrabi, A., & Dominik, M. 2011, MNRAS, 410, 912
'''

# speed of light, c, (au/day)
# orbital radius, a, (au)
# mass of source star, m_star, mass of lens, M (solal masses)
# inclination angle with respect to observer-lens line of sight, φ (radians)
# the time of maximum brightness due to lensing, set at 0 hours
# the Einstein Crossing time, the time taken for the source star to cross the black hole's Einstein Radius, t_e (hours)
# the smallest angular seperation between the source star and the black hole, u_0 (unitless)
# the angular seperation between the source star and the black hole in terms of t, u_t (unitless)
# the universal gravitational constant, G (au^3  *  day^-2  *  solar mass^-1)
# the Schwarzchild Radius, the radius of the even horizon of a black hole, R_S (au)


# required inputs: set of time values (t_vals, 100000 values), t_0 (set to 0), φ, M, m_star, a

def rahvar_main(t, t_0, φ, M, m_star, a):

    # constants
    c = 173.145 # (au/day)
    G = 0.000295913010 # au * 1/solar mass * (au/day)^2

    # equation 2
    R_S = (2 * M * G)/(c ** 2) # au

    # equation 17
    t_E = ((2 * a)/c) * (np.sqrt(M/(m_star + M))) * 1/24

    # equation 19
    u_0 = np.sqrt(φ * (a / ((2 * R_S))))

    # equation 2
    u_t = np.sqrt((u_0 ** 2) + (((t - t_0)/t_E) ** 2))

    # equaiton 14
    Amplitude = (2/np.pi) * ((1 + 1/(u_t ** 2)) * np.arcsin(1/(np.sqrt(1 + 1/(u_t ** 2)))) + 1/u_t)

    return Amplitude

def rahvar_init(t_vals, t_0, φ, M, m_star, a): # call to run rahvar_main

    amp_vals = np.zeros(100000)

    additive = 0

    for element in t_vals:
        amp_vals[additive] = rahvar_main(element, t_0, φ, M, m_star, a)
        additive += 1

    return amp_vals # list of amplitude values for each given time


# example inputs below

'''
t_0 = 0.0 # hours
φ = 1.59989 * (10.0 ** (-6)) # radians
M = 8.5 # solar mass
m_star = 0.35 # solar amss
a = 17 # au

t_vals = np.linspace(-2.5, 2.5, 100000) # hours
'''

amp_vals = rahvar_init(t_vals, t_0, φ, M, m_star, a)

plt.figure(1)
plt.plot(t_vals, amp_vals, label="Rahvar")
plt.xlabel("Time [Hours]")
plt.ylabel("Magnfication")

plt.show()
