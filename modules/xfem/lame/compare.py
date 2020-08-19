import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import pandas as pd

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 24

plt.figure(figsize=(24, 8))

plt.subplot(1, 3, 1)
df = pd.read_csv(
    'no_xfem/output/data_nr_96_nt_96_u_r_0011.csv')
plt.plot(df['x'], df['disp_x'], label='analytical solution')
df = pd.read_csv(
    'xfem_finite_strain_heal/output/data_nr_95_nt_95_u_r_0011.csv')
plt.plot(df['x'], df['disp_x'], label='finite strain with healing-and-re-cut')
plt.xlabel('$r$')
plt.ylabel('$u_r$')
plt.legend(loc='upper center', fontsize='small')
plt.grid()

plt.subplot(1, 3, 2)
df = pd.read_csv(
    'no_xfem/output/data_nr_96_nt_96_stress_rr_0011.csv')
plt.plot(df['x'], df['stress_rr'])
df = pd.read_csv(
    'xfem_finite_strain_heal/output/data_nr_95_nt_95_stress_rr_0011.csv')
plt.plot(df['x'], df['stress_rr'])
plt.xlabel('$r$')
plt.ylabel('$\sigma_{rr}$')
plt.grid()

plt.subplot(1, 3, 3)
df = pd.read_csv(
    'no_xfem/output/data_nr_96_nt_96_stress_tt_0011.csv')
plt.plot(df['x'], df['stress_tt'])
df = pd.read_csv(
    'xfem_finite_strain_heal/output/data_nr_95_nt_95_stress_tt_0011.csv')
plt.plot(df['x'], df['stress_tt'])
plt.xlabel('r')
plt.ylabel('$\sigma_{\\theta\\theta}$')
plt.grid()

plt.tight_layout()

plt.show()
