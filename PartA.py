import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
import CubicEquationSolver


R = 8.314e-5  # universal gas constant, m3-bar/K-mol
# for oxygen:
Pc = 50.0
Tc = 154.58
omega = 0.022

# for Nitric Oxide:
Pc_NO = 64.8
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
    A = a*P/(R**2 * T**2)
    B = (b*P)/(R*T)
    liquid_root = CubicEquationSolver.solve(1, (B-1), (A-3*B**2-2*B), (A*B-B**2-B**3))*R*T/P
    return liquid_root[np.where(liquid_root > 0, liquid_root, np.inf).argmax() - 1]


def data_gen(T, Vm):
    Tr_NO = T / Tc_NO
    a = 0.45724 * ((R ** 2 * Tc_NO ** 2) / Pc_NO)
    b = 0.07780 * ((R * Tc_NO) / Pc_NO)
    kappa = 0.37464 + (1.54226 * omega_NO) - (0.26992 * omega_NO ** 2)
    alpha = (1 + kappa * (1 - Tr_NO ** 0.5)) ** 2
    P = ((R*T)/(Vm - b)) - ((a*alpha)/(Vm**2 - 2*b*Vm-b**2))
    return P


def Nitric_Oxide_Graph(T):
    n = 100
    Vm_list = []
    p_list = []
    Vm_list.append(10.0)
    for i in range(n-1):
        Vm_list.append(Vm_list[i-1]*1.2)
    for j in range(n):
        p_list.append(data_gen(T,Vm_list[j]))

    return p_list, Vm_list


def plot_NO_100():
    fig = plt.figure()
    ax = fig.add_subplot(2,1,1)

    ax.plot(Nitric_Oxide_Graph(100)[0], Nitric_Oxide_Graph(100)[1], color='blue', lw=2, label="T=100k")
    ax.set_xscale('log')
    ax.set_xlabel('Molar Volume(m^3/mol)')
    ax.set_ylabel('Pressure(bar)')
    ax.grid(linestyle='--')
    plt.title('Vm vs P for NO')
    plt.legend()

    plt.show()

while True:
    choice = input("Part A"+'\n' + '1.Calculate Pressure for O2 with Temperature and Molar Volume'+ '\n' +
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
        print('The Molar Volume is:', Vm_Calc(temp, P),'m^3/mol')
    elif choice == 3:
        plot_NO_100()

