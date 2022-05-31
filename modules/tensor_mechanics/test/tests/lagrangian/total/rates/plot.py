#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    G = 5000  # shear modulus

    def analytical_jaumann(k):
        return np.sin(k)

    def analytical_green_naghdi(k):
        b = np.arctan(k / 2)
        return 2 * np.cos(2 * b) * (2 * b - 2 * np.tan(2 * b) * np.log(np.cos(b)) - np.tan(b))

    truesdell = np.loadtxt("truesdell_shear_out.csv", skiprows=1,
                           delimiter=',')
    jaumann = np.loadtxt("jaumann_shear_out.csv", skiprows=1,
                         delimiter=',')
    green_naghdi = np.loadtxt("green_naghdi_shear_out.csv", skiprows=1,
                              delimiter=',')

    plt.plot(jaumann[:, 0] / 2, analytical_jaumann(jaumann[:, 0] / 2),
             'r-', label="analytical - Jaumann")
    plt.plot(green_naghdi[:, 0] / 2, analytical_green_naghdi(jaumann[:, 0] / 2),
             'b-', label="analytical - Green-Naghdi")
    plt.plot(truesdell[:, 0] / 2, truesdell[:, 2],
             'ko', markevery=5, label="numerical - Truesdell")
    plt.plot(jaumann[:, 0] / 2, jaumann[:, 2],
             'ro', markevery=5, label="numerical - Jaumann")
    plt.plot(green_naghdi[:, 0] / 2, green_naghdi[:, 2],
             'bo', markevery=5, label="numerical - Green-Naghdi")
    plt.xlabel("Shear strain")
    plt.ylabel("Shear stress")
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()
