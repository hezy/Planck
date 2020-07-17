"""
# -*- coding: utf-8 -*-
Created on Apr 30, 2019
@author: Hezy Amiel
planck.py
this script uses Planck's equation to create a set of black body radiation
curves at various temperatures, for each curve find λmax
(the wavelength of maximum radiation)
and use it to plot Wien's law: λmax = b/T
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import constants
import pandas as pd


# Physical constants
h = constants.value(u'Planck constant')
k = constants.value(u'Boltzmann constant')
c = constants.value(u'speed of light in vacuum')
b = constants.value(u'Wien wavelength displacement law constant')


def planck(lam, T):
    '''
    Palnck's function -
    T = temperature (K), lam = wavelength (nm), l = wavelength (m)
    this calculates the spectral radiance - the power per unit solid angle and
    per unit of area normal to the propagation, density of wavelength λ
    radiation per unit wavelegth at thermal equilibrium at temperature T.
    '''
    # converting nanometer to meter
    l = lam * 1e-9
    # Planck's function
    return (8*np.pi*h*c/l**5)/(np.exp(h*c/(l*k*T)-1))


def noise(signal, sigma_background, sigma_measurment):
    '''
    Adds noise to a given signal
    sigma_background is the standart deviation of a normal background noise
    sigam_measurment is the standart deviation of a normal noise proportional
    to the signal
    '''
    noisy_signal = (signal
                    + np.random.normal(0, sigma_background, np.size(signal))
                    + signal * np.random.normal(0,
                                                sigma_measurment,
                                                np.size(signal)))
    return noisy_signal


def wien(x, a, b):
    '''
    a function for fitting Wien's curve λmax = b/T
    a is due to approximation error, b = Wien's constant (nm*K),
    to be evaluated. expected value is b = 2.8977729e−3 m*K
    '''
    return a + b/x


def write_file(file_name, title, x, y):
    with open(file_name, 'w') as file:
        file.write(title+'\n')
    df = pd.DataFrame({'wavelength (nm)': x,
                       'Intesity (arb)': y})
    df.to_csv(file_name,
              mode='a',
              index=False,
              header=True,
              sep='\t')


# preconfiguring the figure
plt.close('all')
plt.rcParams.update({'font.size': 12})
fig, axs = plt.subplots(nrows=2, ncols=1, sharex=False, figsize=(8, 12))
plt.subplots_adjust(left=None, bottom=None, right=None,
                    top=0.95, wspace=None, hspace=0.3)
fst = 16  # font size for title
fsl = 14  # font size for axes labels


# create Planck's curves for different temperatures
x = np.arange(400, 2600, 100)  # wavelength (nm)
xfit = np.arange(1, 10000, 1)  # wavelength (nm)
T = np.arange(1500, 5500, 500)  # temperature (K)
Lmax = np.array([])

R0 = 1.12
Rw = 0.32
alpha = 4.5e-3

for t in T:
    # R
    R = R0 * (1 + t/293) + Rw

    # creating Planck's curve for temperature = t (with some random noise)
    y = (planck(x, t))  # * (1 + np.random.normal(0, 0.1, None))
    y = noise(y, 0, 0.02)
    axs[0].plot(x, y, "o")

    # fitting Planck's curve for t
    popt, pcov = curve_fit(planck, x, y, 2000)
    print('fitted T =', popt, '+/-', pcov)
    yfit = planck(xfit, *popt)
    axs[0].plot(xfit, yfit, "-")

    #write_file(str(t) + 'K.csv', 'T = ' + str(t) + 'K', x, y)
    write_file(str(R) + 'Ohm.csv', 'T = ' + str(R) + 'Ohm', x, y)

    '''
    finding the wavelength of maximal radiance for Planck's curve and
    appending to Wien's curve
    '''
    L = np.argmax(yfit)
    Lmax = np.append(Lmax, np.argmax(yfit))  # TODO


# fitting Wien's curve with a + b/T
axs[1].plot(T, Lmax, "bo")
Tfit = np.arange(100, 6000, 100)  # temperature (K)
popt, pcov = curve_fit(wien, T, Lmax)
print('a,b = ')
print(popt)
print('Delta a,b = ')
print(pcov)
axs[1].plot(Tfit, wien(Tfit, *popt), "r-")


# aranging graph 1
axs[0].grid(True)
axs[0].set_title("Planck's curves", fontsize=fst)
axs[0].set_xlabel("wavelength (nm)", fontsize=fsl)
axs[0].set_xlim(200, 2500)
axs[0].set_ylabel(r'$spectral \: radiance \: (W/sr/m^3)$', fontsize=fsl)
# aranging graph 2
axs[1].grid(True)
axs[1].set_title("Wien’s Law", fontsize=fst)
axs[1].set_xlabel("Temperature (K)", fontsize=fsl)
axs[1].set_xlim(1, 6000)
axs[1].set_ylabel(r'$\lambda_{max} \: (nm)$', fontsize=fsl)
axs[1].set_ylim(0, 10000)
# show figure
plt.show()
