import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import pandas as pd

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 24

plt.figure(figsize=(18, 6))

plt.subplot(1, 3, 1)
df = pd.read_csv(
    'stress.csv')
plt.plot(df['Points:0'], df['stress_rr'], label='MOOSE')
df = pd.read_csv(
    'benchmark_stress_rr.csv', header=None)
plt.plot(df[0], df[1], label='analytical')
plt.xlim(13.8, 25)
plt.xlabel('$r$')
plt.ylabel('$\sigma_{rr}$')
plt.legend(loc='upper center', fontsize='small')
plt.grid()

plt.subplot(1, 3, 2)
df = pd.read_csv(
    'stress.csv')
plt.plot(df['Points:0'], df['stress_zz'], label='MOOSE')
df = pd.read_csv(
    'benchmark_stress_zz.csv', header=None)
plt.plot(df[0], df[1], label='analytical')
plt.xlim(13.8, 16)
plt.xlabel('$r$')
plt.ylabel('$\sigma_{zz}$')
plt.grid()

plt.subplot(1, 3, 3)
df = pd.read_csv(
    'stress.csv')
plt.plot(df['Points:0'], df['stress_tt'], label='MOOSE')
df = pd.read_csv(
    'benchmark_stress_tt.csv', header=None)
plt.plot(df[0], df[1])
plt.xlim(13.8, 16)
plt.xlabel('$r$')
plt.ylabel('$\sigma_{\\theta\\theta}$')
plt.grid()

plt.tight_layout()

plt.figure(figsize=(6, 6))

A = []
stress_rr_max = []
for a in np.linspace(0, 0.15, 16):
    df = pd.read_csv(
        'stress_wavy_A_' + "{:.2f}".format(a) + '.csv')
    A.append(a * 1000)
    stress_rr_max.append(df['stress_rr_max'][0])
plt.plot(A, stress_rr_max, 'o-', label='MOOSE')
df = pd.read_csv(
    'benchmark_A.csv', header=None)
plt.plot(df[0], df[1], label='COMSOL')
plt.xlim(0, 150)
plt.xlabel('$A$')
plt.ylabel('$\sigma_{rr,\mathrm{max}}$')
plt.legend(loc='lower right', fontsize='small')
plt.grid()

plt.tight_layout()

plt.show()
