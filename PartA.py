import sys

import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt
from shapely.geometry import LineString
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib.ticker import FormatStrFormatter, FixedLocator
import CubicEquationSolver

n = 1000000
R = 8.314  # universal gas constant, m3-bar/K-mol
# for oxygen:
Pc = 50.0
Tc = 154.58
omega = 0.022

# for Nitric Oxide:
Pc_NO = 6.485e6
Tc_NO = 180.2
omega_NO = 0.583


def pressure_Calc(T, Vm):
    Tr = T / Tc
    a = 0.45724 * ((R ** 2 * Tc ** 2) / Pc)
    b = 0.07780 * ((R * Tc) / Pc)
    kappa = 0.37464 + (1.54226 * omega) - (0.26992 * omega ** 2)
    alpha = (1 + kappa * (1 - Tr ** 0.5)) ** 2

    p = ((R * T) / (Vm - b)) - (a * alpha / (Vm ** 2 + 2 * b * Vm - b ** 2))
    return p


def Vm_Calc(T, P):
    a = 0.45724 * ((R ** 2 * Tc ** 2) / Pc)
    b = 0.07780 * ((R * Tc) / Pc)
    A = a * P / (R ** 2 * T ** 2)
    B = (b * P) / (R * T)
    print(A, B)
    liquid_root = CubicEquationSolver.solve(1, (B - 1), (A - 3 * B ** 2 - 2 * B), (A * B - B ** 2 - B ** 3)) * R * T / P
    return liquid_root  # [np.where(liquid_root > 0, liquid_root, np.inf).argmax() - 1]


# Solves for the pressures of an T degree Isotherm given a set of molar volumes
def data_gen(T, Vm):
    Tr_NO = T / Tc_NO
    a = 0.45724 * ((R ** 2 * Tc_NO ** 2) / Pc_NO)
    b = 0.07780 * ((R * Tc_NO) / Pc_NO)
    kappa = 0.37464 + (1.54226 * omega_NO) - (0.26992 * omega_NO ** 2)
    alpha = (1 + kappa * (1 - Tr_NO ** 0.5)) ** 2
    P = ((R * T) / (Vm - b)) - ((a * alpha) / (Vm ** 2 + 2 * b * Vm - b ** 2))
    return P


def data_gen2(Vm):
    T = 100
    Tr_NO = T / Tc_NO
    a = 0.45724 * ((R ** 2 * Tc_NO ** 2) / Pc_NO)
    b = 0.07780 * ((R * Tc_NO) / Pc_NO)
    kappa = 0.37464 + (1.54226 * omega_NO) - (0.26992 * omega_NO ** 2)
    alpha = (1 + kappa * (1 - Tr_NO ** 0.5)) ** 2
    P = ((R * T) / (Vm - b)) - ((a * alpha) / (Vm ** 2 + 2 * b * Vm - b ** 2))
    return P


# Obtains values to plot NO graph, including P_min and P_max
def Nitric_Oxide_Graph(T):
    # Initializing function molar volume (x) and pressure values (y)

    Vm_list = np.linspace(2.0e-5, 50, n)
    p_list = np.zeros(n)


    # Solves for pressure values using data_gen
    for i in range(n):
        p_list[i] = data_gen(T, Vm_list[i])
    return Vm_list, p_list

def Nitric_Oxide_Graph_Visual(T):
    # Initializing function molar volume (x) and pressure values (y)

    Vm_list = np.linspace(2.0e-5, 1, n)
    p_list = np.zeros(n)


    # Solves for pressure values using data_gen
    for i in range(n):
        p_list[i] = data_gen(T, Vm_list[i])
    return Vm_list, p_list


def plot_NO_100(visual):
    if visual:
        x, y = Nitric_Oxide_Graph_Visual(100)
    else:
        x, y = Nitric_Oxide_Graph(100)
    fig = plt.figure()
    fig.set_size_inches(10.5, 10.5)
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(x, y, color='blue', lw=2, label="T=100k")
    # idx = np.argwhere(np.diff(np.sign(g-y))).flatten()
    # ax.plot(x[idx], y[idx], 'ro')
    ax.set_xlabel('Molar Volume(m^3/mol)')
    ax.set_ylabel('Pressure(Pa)')
    ax.set_xscale('log')
    ax.grid(linestyle='--')
    plt.title('Vm vs P for NO')
    plt.legend()

    plt.show()




def big_chungus(pressure_guess, x, y):
    f = np.linspace(1.0e-5, 1.0, 1000000)
    g = np.full(len(f), pressure_guess)
    line_1 = LineString(np.column_stack((f, g)))
    line_2 = LineString(np.column_stack((x, y)))
    intersection = line_1.intersection(line_2)
    print(intersection)
    mom = integrate.quad(data_gen2, intersection[0].x, intersection[1].x)
    dad = integrate.quad(data_gen2, intersection[1].x, intersection[2].x)

    return abs(mom[0]), abs(dad[0])


def Guess_n_check(pressure_guesses, x, y):
    i = int(len(pressure_guesses)/2)
    AreaI, AreaII = big_chungus(pressure_guesses[i],x, y)
    if abs(1 - AreaI / AreaII) < 0.01:
        return pressure_guesses[i]
    while True:
        if AreaI > AreaII:
            i -= 1
            AreaI, AreaII = big_chungus(pressure_guesses[i],x,y)
            if abs(1 - AreaI / AreaII) < 0.01:
                return pressure_guesses[i]
            elif AreaII > AreaI:
                new_pressures = np.linspace(pressure_guesses[i], pressure_guesses[i + 1], 100)
                print('ahh')
                Guess_n_check(new_pressures,x,y)
        elif AreaII > AreaI:
            i += 1
            AreaI, AreaII = big_chungus(pressure_guesses[i],x,y)
            if abs(1 - AreaI / AreaII) < 0.01:
                return pressure_guesses[i]
            elif AreaI > AreaII:
                new_pressures = np.linspace(pressure_guesses[i - 1], pressure_guesses[i], 100)
                print('ahh')
                Guess_n_check(new_pressures,x,y)


while True:
    choice = input("Part A" + '\n' + '1.Calculate Pressure for O2 with Temperature and Molar Volume' + '\n' +
                   '2. Calculate molar volume for O2 using Temperature and Pressure' + '\n' + '3. Display the Vm vs P plot for Nitric Oxide at T = 100k')
    choice = int(choice)
    if choice == 1:
        temp = input('Enter the temperature:')
        temp = float(temp)
        Vm = input('Enter the Molar Volume:')
        Vm = float(Vm)
        print('The Pressure is:', pressure_Calc(temp, Vm), 'bar')
    elif choice == 2:
        temp = input('Enter the temperature:')
        temp = float(temp)
        P = input('Enter the Pressure:')
        P = float(P)
        print('The Molar Volume is:', Vm_Calc(temp, P), 'm^3/mol')
    elif choice == 3:
        plot_NO_100(True)

    elif choice == 4:
        x,y = Nitric_Oxide_Graph(100)
        guesses = np.linspace(71000, 729231, 100)
        ahaha = Guess_n_check(guesses,x,y)
        x, y = Nitric_Oxide_Graph_Visual(100)
        fig = plt.figure()
        fig.set_size_inches(10.5, 10.5)
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(x, y, color='blue', lw=2, label="T=100k")
        # idx = np.argwhere(np.diff(np.sign(g-y))).flatten()
        # ax.plot(x[idx], y[idx], 'ro')
        ax.set_xlabel('Molar Volume(m^3/mol)')
        ax.set_ylabel('Pressure(Pa)')
        ax.set_xscale('log')
        ax.grid(linestyle='--')
        ax.set_xlim(1e-5, 0.004)
        plt.title('Vm vs P for NO')
        plt.legend()
        f = np.linspace(1.0e-5, 50, 1000000)
        g = np.full(len(f), ahaha)
        ax.plot(f,g, color='red', lw=2)
        plt.show()




    elif choice == 5:
        x, y, a, b = Nitric_Oxide_Graph(100)
        z = input('whan is big chingus')
        print(big_chungus(z,x,y))

