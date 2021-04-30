import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter, FixedLocator
import CubicEquationSolver

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
    Tr = T / Tc_NO
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
    liquid_root = CubicEquationSolver.solve(1, (B - 1), (A - 3 * B ** 2 - 2 * B), (A * B - B ** 2 - B ** 3)) * R * T / P
    return liquid_root[np.where(liquid_root > 0, liquid_root, np.inf).argmax() - 1]


#Solves for the pressures of an T degree Isotherm given a set of molar volumes
def data_gen(T, Vm):
    Tr_NO = T / Tc_NO
    a = 0.45724 * ((R ** 2 * Tc_NO ** 2) / Pc_NO)
    b = 0.07780 * ((R * Tc_NO) / Pc_NO)
    kappa = 0.37464 + (1.54226 * omega_NO) - (0.26992 * omega_NO ** 2)
    alpha = (1 + kappa * (1 - Tr_NO ** 0.5)) ** 2
    P = ((R * T) / (Vm - b)) - ((a * alpha) / (Vm ** 2 + 2 * b * Vm - b ** 2))
    return P

#Obtains values to plot NO graph, including P_min and P_max
def Nitric_Oxide_Graph(T):
    #Initializing function molar volume (x) and pressure values (y)
    n = 1000000
    Vm_list = np.linspace(2.0e-5, 1.0,n)
    p_list = np.zeros(n)

    #Initializing min/max values and indices
    P_min = 0
    P_max = 0
    P_minIndex = 0
    P_maxIndex = 0
    
    #Solves for pressure values using data_gen
    for i in range (n):
        p_list[i] = data_gen(T, Vm_list[i])
        if(p_list[i-1] < 0 and p_list[i] > 0):
            print(p_list[i-1], " ",p_list[i])
        if(p_list[i] < P_min):
            P_min = p_list[i]
            P_minIndex = i
            print(P_minIndex)
    
    #Sorting to find maximum
    for i in range (P_minIndex,n):

        if(p_list[i] > P_max):
            P_max = p_list[i]
            P_maxIndex = i

    #Creates coordinate points for min and max
    a = [Vm_list[P_minIndex], P_min]
    b = [Vm_list[P_maxIndex], P_max]

    
        
    return Vm_list, p_list,a,b
    

def plot_NO_100():
    x,y,a,b  = Nitric_Oxide_Graph(100)
    fig = plt.figure()
    fig.set_size_inches(10.5, 10.5)
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(x, y, color='blue', lw=2, label="T=100k")
    ax.set_xlabel('Molar Volume(m^3/mol)')
    ax.set_ylabel('Pressure(Pa)')
    ax.set_xscale('log')
    ax.grid(linestyle='--')
    plt.title('Vm vs P for NO')
    plt.legend()

    plt.draw()

    return a,b


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
        a,b = plot_NO_100()
        print(a[0], " ", a[1])
        print(b[0], " ", b[1])

        #Part B) Using the algorithm to find P_vap
        #Beginning with a guess of 0
        p_guess = 0.0
        

        plt.show()

        