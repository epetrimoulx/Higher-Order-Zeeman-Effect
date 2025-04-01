import matplotlib.pyplot as plt
import numpy as np 
import seaborn as sns

FINE_STRUCTURE_CONST = 1.0/137.0
H_GROUND_STATE_ENERGY = 0.5
PROTON_MASS_AU = 1836.15267
NUCLEAR_MAGNETON = 1 / (2 * PROTON_MASS_AU)
NUCLEAR_G_FACTOR = -4.25509
BOHR_MAGNETON = 1/2
ELECTRONIC_G_FACTOR = 2
NUCLEAR_CHARGE = 2

def Calc_Non_Relativistic_Energy():
    return -H_GROUND_STATE_ENERGY 
    
def Calc_Nuclear_Zeeman_Energy(Magnetic_Field, Spin_Magnetic_Quantum_Number_Nucleus):
    return NUCLEAR_MAGNETON * Magnetic_Field * NUCLEAR_G_FACTOR * Spin_Magnetic_Quantum_Number_Nucleus

def Calc_Electronic_Zeeman_Energy(Magnetic_Field, Spin_Magnetic_Quantum_Number_electron):
    return BOHR_MAGNETON * Magnetic_Field * ELECTRONIC_G_FACTOR * Spin_Magnetic_Quantum_Number_electron

def Calc_Relativistic_Energy(Magnetic_Field, Spin_Magnetic_Quantum_Number_electron):
    return -0.5 * Magnetic_Field**3 * NUCLEAR_CHARGE**(-7/2) * FINE_STRUCTURE_CONST**2 * Spin_Magnetic_Quantum_Number_electron

def Graph_Energy(Magnetic_Field, Total_Energy, Labels) -> None:
    sns.set_theme(style="whitegrid") 
    colors = sns.color_palette("Set2", n_colors=4) 
    linestyles = ['-', '--', '-.', ':']
    plt.figure(figsize=(12,8))
    plt.title('Higher-order Zeeman Effect Corrections', fontsize=18, weight='bold')
    plt.xlabel(r'$B$ Field [T]', fontsize=14)
    plt.ylabel('Energy [a.u.]', fontsize=14)

    for i, (energy, label, color, ls) in enumerate(zip(Total_Energy, Labels, colors, linestyles)):
        sns.lineplot(x=Magnetic_Field * 2.3505175e5, y=energy, label=label, color=color, linestyle=ls, linewidth=2)
    
    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
    plt.legend(loc='upper left', fontsize=12, bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.tick_params(axis='both', direction='in', length=6, width=1.5, colors='black', grid_color='gray', grid_alpha=0.7)
    plt.show()
    
def Graph_Ratio(Magnetic_Field, Energy_Ratio, Labels) -> None:
    # GRAPH 1 - NO LOG SCALE
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(12, 8))
    plt.title(r'$\frac{C_{rel}^{(2)}}{E_{nuclear}}$ Ratio', fontsize=18, weight='bold', y=1.02)
    plt.xlabel(r'$B$ Field [T]', fontsize=14)
    plt.ylabel(r'$\frac{C_{rel}^{(2)}}{E_{nuclear}}$', fontsize=14, rotation=0, labelpad=20)
    
    linestyles = ['-', '--', '-.', ':']
    colors = sns.color_palette("Set2", n_colors=4)
    
    for i, (energy, label, color, ls) in enumerate(zip(Energy_Ratio, Labels, colors, linestyles)):
        sns.lineplot(x=Magnetic_Field * 2.3505175e5, y=energy, label=label, color=color, linestyle=ls, linewidth=2)

    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
    plt.legend(loc='upper left', fontsize=12, bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.show()
    
    # GRAPH 2 BOTH AXIS LOG SCALE
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(12, 8))
    plt.title(r'$\frac{C_{rel}^{(2)}}{E_{nuclear}}$ Ratio', fontsize=18, weight='bold', y=1.02)
    plt.xlabel(r'$B$ Field [T]', fontsize=14)
    plt.ylabel(r'$\frac{C_{rel}^{(2)}}{E_{nuclear}}$', fontsize=14, rotation=0, labelpad=20)
    plt.xscale('linear')
    plt.yscale('log')
    
    linestyles = ['-', '--', '-.', ':']
    colors = sns.color_palette("Set2", n_colors=4)
    
    for i, (energy, label, color, ls) in enumerate(zip(Energy_Ratio, Labels, colors, linestyles)):
        sns.lineplot(x=Magnetic_Field*2.3505175e5, y=energy, label=label, color=color, linestyle=ls, linewidth=2)

    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
    plt.legend(loc='upper left', fontsize=12, bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.show()
    
def main() -> None:
    Magnetic_Field = np.linspace(1, 50, 1024)
    Magnetic_Field = Magnetic_Field / (2.3505175e5)

    Labels = [
        r'$m_I=-\frac{1}{2}, m_s=-\frac{1}{2}$',
        r'$m_I=-\frac{1}{2}, m_s=+\frac{1}{2}$',
        r'$m_I=+\frac{1}{2}, m_s=-\frac{1}{2}$',
        r'$m_I=+\frac{1}{2}, m_s=+\frac{1}{2}$'
    ]
    
    Total_Energy = [Calc_Non_Relativistic_Energy() + Calc_Nuclear_Zeeman_Energy(Magnetic_Field, -0.5) + Calc_Electronic_Zeeman_Energy(Magnetic_Field, -0.5) + Calc_Relativistic_Energy(Magnetic_Field, -0.5) - 4e-6,
                    Calc_Non_Relativistic_Energy() + Calc_Nuclear_Zeeman_Energy(Magnetic_Field, -0.5) + Calc_Electronic_Zeeman_Energy(Magnetic_Field, 0.5) + Calc_Relativistic_Energy(Magnetic_Field, 0.5),
                    Calc_Non_Relativistic_Energy() + Calc_Nuclear_Zeeman_Energy(Magnetic_Field, 0.5) + Calc_Electronic_Zeeman_Energy(Magnetic_Field, -0.5) + Calc_Relativistic_Energy(Magnetic_Field, -0.5),
                    Calc_Non_Relativistic_Energy() + Calc_Nuclear_Zeeman_Energy(Magnetic_Field, 0.5) + Calc_Electronic_Zeeman_Energy(Magnetic_Field, 0.5) + Calc_Relativistic_Energy(Magnetic_Field, 0.5) + 4e-6
    ]
    
    Graph_Energy(Magnetic_Field, Total_Energy, Labels)
    
    Energy_Ratio = [Calc_Relativistic_Energy(Magnetic_Field, -0.5) / Calc_Nuclear_Zeeman_Energy(Magnetic_Field, -0.5),
                    Calc_Relativistic_Energy(Magnetic_Field, 0.5) / Calc_Nuclear_Zeeman_Energy(Magnetic_Field, -0.5),
                    Calc_Relativistic_Energy(Magnetic_Field, -0.5) / Calc_Nuclear_Zeeman_Energy(Magnetic_Field, 0.5),
                    Calc_Relativistic_Energy(Magnetic_Field, 0.5) / Calc_Nuclear_Zeeman_Energy(Magnetic_Field, 0.5) 
    ]
    
    Graph_Ratio(Magnetic_Field, Energy_Ratio, Labels)

    Energy_Ratio_455T = [Calc_Relativistic_Energy(45.5/2.3505175e5, -0.5) / Calc_Nuclear_Zeeman_Energy(45.5/2.3505175e5, -0.5),
                    Calc_Relativistic_Energy(45.5/2.3505175e5, 0.5) / Calc_Nuclear_Zeeman_Energy(45.5/2.3505175e5, -0.5),
                    Calc_Relativistic_Energy(45.5/2.3505175e5, -0.5) / Calc_Nuclear_Zeeman_Energy(45.5/2.3505175e5, 0.5),
                    Calc_Relativistic_Energy(45.5/2.3505175e5, 0.5) / Calc_Nuclear_Zeeman_Energy(45.5/2.3505175e5, 0.5)
                    ]

    Energy_Ratio_145T = [
        Calc_Relativistic_Energy(14.5 / 2.3505175e5, -0.5) / Calc_Nuclear_Zeeman_Energy(14.5 / 2.3505175e5, -0.5),
        Calc_Relativistic_Energy(14.5 / 2.3505175e5, 0.5) / Calc_Nuclear_Zeeman_Energy(14.5 / 2.3505175e5, -0.5),
        Calc_Relativistic_Energy(14.5 / 2.3505175e5, -0.5) / Calc_Nuclear_Zeeman_Energy(14.5 / 2.3505175e5, 0.5),
        Calc_Relativistic_Energy(14.5 / 2.3505175e5, 0.5) / Calc_Nuclear_Zeeman_Energy(14.5 / 2.3505175e5, 0.5)
        ]

    Energy_Ratio_570T = [
        Calc_Relativistic_Energy(5.7 / 2.3505175e5, -0.5) / Calc_Nuclear_Zeeman_Energy(5.7 / 2.3505175e5, -0.5),
        Calc_Relativistic_Energy(5.7 / 2.3505175e5, 0.5) / Calc_Nuclear_Zeeman_Energy(5.7 / 2.3505175e5, -0.5),
        Calc_Relativistic_Energy(5.7 / 2.3505175e5, -0.5) / Calc_Nuclear_Zeeman_Energy(5.7 / 2.3505175e5, 0.5),
        Calc_Relativistic_Energy(5.7 / 2.3505175e5, 0.5) / Calc_Nuclear_Zeeman_Energy(5.7 / 2.3505175e5, 0.5)
    ]

    Energy_Ratio_117T = [
        Calc_Relativistic_Energy(11.7 / 2.3505175e5, -0.5) / Calc_Nuclear_Zeeman_Energy(11.7 / 2.3505175e5, -0.5),
        Calc_Relativistic_Energy(11.7 / 2.3505175e5, 0.5) / Calc_Nuclear_Zeeman_Energy(11.7 / 2.3505175e5, -0.5),
        Calc_Relativistic_Energy(11.7 / 2.3505175e5, -0.5) / Calc_Nuclear_Zeeman_Energy(11.7 / 2.3505175e5, 0.5),
        Calc_Relativistic_Energy(11.7 / 2.3505175e5, 0.5) / Calc_Nuclear_Zeeman_Energy(11.7 / 2.3505175e5, 0.5)
    ]

    print(Energy_Ratio_455T)
    print(Energy_Ratio_145T)
    print(Energy_Ratio_570T)
    print(Energy_Ratio_117T)

    Magnetic_Field = 45.5 / 2.3505175e5
    Total_Energy_New = [
        Calc_Non_Relativistic_Energy() + Calc_Nuclear_Zeeman_Energy(Magnetic_Field,-0.5) + Calc_Electronic_Zeeman_Energy(Magnetic_Field, -0.5) + Calc_Relativistic_Energy(Magnetic_Field, -0.5),
        Calc_Non_Relativistic_Energy() + Calc_Nuclear_Zeeman_Energy(Magnetic_Field,-0.5) + Calc_Electronic_Zeeman_Energy(Magnetic_Field, 0.5) + Calc_Relativistic_Energy(Magnetic_Field, 0.5),
        Calc_Non_Relativistic_Energy() + Calc_Nuclear_Zeeman_Energy(Magnetic_Field,  0.5) + Calc_Electronic_Zeeman_Energy(Magnetic_Field, -0.5) + Calc_Relativistic_Energy(Magnetic_Field, -0.5),
        Calc_Non_Relativistic_Energy() + Calc_Nuclear_Zeeman_Energy(Magnetic_Field,0.5) + Calc_Electronic_Zeeman_Energy(Magnetic_Field, 0.5) + Calc_Relativistic_Energy(Magnetic_Field, 0.5)
    ]

    Total_Energy_Old = [
        Calc_Non_Relativistic_Energy() + Calc_Nuclear_Zeeman_Energy(Magnetic_Field,-0.5) + Calc_Electronic_Zeeman_Energy(Magnetic_Field, -0.5),
        Calc_Non_Relativistic_Energy() + Calc_Nuclear_Zeeman_Energy(Magnetic_Field,-0.5) + Calc_Electronic_Zeeman_Energy(Magnetic_Field, 0.5),
        Calc_Non_Relativistic_Energy() + Calc_Nuclear_Zeeman_Energy(Magnetic_Field,0.5) + Calc_Electronic_Zeeman_Energy(Magnetic_Field, -0.5),
        Calc_Non_Relativistic_Energy() + Calc_Nuclear_Zeeman_Energy(Magnetic_Field, 0.5) + Calc_Electronic_Zeeman_Energy(Magnetic_Field, 0.5)
    ]

    print(Total_Energy_New)

    for i in range(len(Total_Energy_Old)):
        print(Total_Energy_Old[i] - Total_Energy_New[i])


if __name__ == '__main__':
    main()