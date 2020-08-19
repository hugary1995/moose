import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import pandas as pd

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 24

plt.figure(figsize=(24, 12))

plt.subplot(2, 3, 1)
for n in range(10):
    df = pd.read_csv(
        'xfem_small_strain_moving/output/data_nr_95_nt_95_u_r_' + '{0:04d}'.format(n + 11) + '.csv')
    plt.plot(df['x'], df['disp_x'], label='$r_c$ = ' + str((n + 1) * 0.1 + 2))
plt.xlabel('$r$')
plt.ylabel('$u_r$')
plt.legend(loc='upper right', fontsize='small')
plt.grid()

plt.subplot(2, 3, 2)
for n in range(10):
    df = pd.read_csv(
        'xfem_small_strain_moving/output/data_nr_95_nt_95_stress_rr_' + '{0:04d}'.format(n + 11) + '.csv')
    plt.plot(df['x'], df['stress_rr'])
plt.xlabel('$r$')
plt.ylabel('$\sigma_{rr}$')
plt.title('analytical solution under small strain assumption', pad=20)
plt.grid()

plt.subplot(2, 3, 3)
for n in range(10):
    df = pd.read_csv(
        'xfem_small_strain_moving/output/data_nr_95_nt_95_stress_tt_' + '{0:04d}'.format(n + 11) + '.csv')
    plt.plot(df['x'], df['stress_tt'])
plt.xlabel('r')
plt.ylabel('$\sigma_{\\theta\\theta}$')
plt.grid()

plt.subplot(2, 3, 4)
for n in range(10):
    df = pd.read_csv(
        'xfem_finite_strain_moving/output/data_nr_95_nt_95_u_r_' + '{0:04d}'.format(n + 11) + '.csv')
    plt.plot(df['x'], df['disp_x'])
plt.xlabel('$r$')
plt.ylabel('$u_r$')
plt.grid()

plt.subplot(2, 3, 5)
for n in range(10):
    df = pd.read_csv(
        'xfem_finite_strain_moving/output/data_nr_95_nt_95_stress_rr_' + '{0:04d}'.format(n + 11) + '.csv')
    plt.plot(df['x'], df['stress_rr'])
plt.xlabel('r')
plt.ylabel('$\sigma_{\\theta\\theta}$')
plt.title('finite strain solution with healing-and-re-cut', pad=20)
plt.grid()

plt.subplot(2, 3, 6)
for n in range(10):
    df = pd.read_csv(
        'xfem_finite_strain_moving/output/data_nr_95_nt_95_stress_tt_' + '{0:04d}'.format(n + 11) + '.csv')
    plt.plot(df['x'], df['stress_tt'])
plt.xlabel('r')
plt.ylabel('$\sigma_{\\theta\\theta}$')
plt.grid()

plt.tight_layout()
plt.show()
