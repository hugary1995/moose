import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import pandas as pd

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 24

plt.subplot(1, 1, 1)

df = pd.read_csv(
    'convergence.csv', delimiter=' ')
nx = pd.to_numeric(df['nx'], downcast='float')
ur = pd.to_numeric(df['disp_r'], downcast='float')
plt.plot(nx, ur, '.-', label='no XFEM')

df = pd.read_csv(
    'convergence_xfem.csv', delimiter=' ')
nx = pd.to_numeric(df['nx'], downcast='float')
ur = pd.to_numeric(df['disp_r'], downcast='float')
plt.plot(nx, ur, '.-', label='XFEM without healing-and-re-cut')

df = pd.read_csv(
    'convergence_xfem_heal.csv', delimiter=' ')
nx = pd.to_numeric(df['nx'], downcast='float')
ur = pd.to_numeric(df['disp_r'], downcast='float')
plt.plot(nx, ur, '.-', label='XFEM with healing-and-re-cut')

plt.xlabel('number of elements')
plt.ylabel('$u_r$')
plt.legend(loc='upper center', fontsize='small')
plt.grid()

plt.tight_layout()

plt.show()
