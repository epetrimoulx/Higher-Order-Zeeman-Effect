import mpmath as mp
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

mp.dps = 100  # Set decimal places to 25

FINE_STRUCTURE_CONST = mp.mpf(1) / mp.mpf(137)
H_GROUND_STATE_ENERGY = mp.mpf(0.5)
PROTON_MASS_AU = mp.mpf('1836.15267')
NUCLEAR_MAGNETON = mp.mpf(1) / (2 * PROTON_MASS_AU)
NUCLEAR_G_FACTOR = mp.mpf('-4.25509')
BOHR_MAGNETON = mp.mpf(1) / 2
ELECTRONIC_G_FACTOR = mp.mpf(2)
NUCLEAR_CHARGE = mp.mpf(2)


def Graph_Energy(Magnetic_Field, Total_Energy, Labels) -> None:
    sns.set_theme(style="whitegrid")
    colors = sns.color_palette("Set2", n_colors=4)
    linestyles = ['-', '--', '-.', ':']
    plt.figure(figsize=(12, 8))
    plt.title('Higher-order Zeeman Effect Corrections', fontsize=18, weight='bold')
    plt.xlabel(r'$B$ Field [T]', fontsize=14)
    plt.ylabel('Energy [a.u.]', fontsize=14)
    plt.xlim(49.999999999, 49.9999999999)
    plt.ylim(0.499, 0.499999)

    for i, (energy, label, color, ls) in enumerate(zip(Total_Energy, Labels, colors, linestyles)):
        sns.lineplot(x=list(map(float, Magnetic_Field)), y=list(map(float, energy)), label=label, color=color,
                     linestyle=ls, linewidth=2)

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
        sns.lineplot(x=list(map(float, Magnetic_Field)), y=list(map(float, energy)), label=label, color=color,
                     linestyle=ls, linewidth=2)

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
        sns.lineplot(x=list(map(float, Magnetic_Field)), y=list(map(float, energy)), label=label, color=color,
                     linestyle=ls, linewidth=2)

    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
    plt.legend(loc='upper left', fontsize=12, bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.show()

def ratio(Magnetic_Field, electron_spin, nuclear_spin):
    return (FINE_STRUCTURE_CONST**2 * Magnetic_Field**2 * PROTON_MASS_AU * electron_spin) / (NUCLEAR_G_FACTOR * nuclear_spin * NUCLEAR_CHARGE**mp.mpf('-7/2'))

def ordinary(Magnetic_Field, electron_spin, nuclear_spin):
    return (H_GROUND_STATE_ENERGY + (mp.mpf(1) / (2 * PROTON_MASS_AU)) * NUCLEAR_G_FACTOR * nuclear_spin * Magnetic_Field
            + mp.mpf('0.5') * ELECTRONIC_G_FACTOR * electron_spin * Magnetic_Field)

def New(Magnetic_Field, electron_spin, nuclear_spin):
    relativistic_contribution = (mp.mpf('0.5') * FINE_STRUCTURE_CONST**2 * Magnetic_Field**3 * NUCLEAR_CHARGE**mp.mpf('-7/2') * electron_spin)
    return ordinary(Magnetic_Field, electron_spin, nuclear_spin) + relativistic_contribution

def main():
    Labels = [
        r'$m_I=-\frac{1}{2}, m_s=-\frac{1}{2}$',
        r'$m_I=-\frac{1}{2}, m_s=+\frac{1}{2}$',
        r'$m_I=+\frac{1}{2}, m_s=-\frac{1}{2}$',
        r'$m_I=+\frac{1}{2}, m_s=+\frac{1}{2}$'
    ]

    fields = [1.45, 5.7, 11.7, 45.5]
    for B in fields:
        B_scaled = mp.mpf(B) / mp.mpf('2.3505175e5')
        print(f"\n{B} T")
        print("\nRATIOS")
        print(mp.nstr(ratio(B_scaled, 0.5, 0.5), 25))
        print(mp.nstr(ratio(B_scaled, 0.5, -0.5), 25))
        print(mp.nstr(ratio(B_scaled, -0.5, 0.5), 25))
        print(mp.nstr(ratio(B_scaled, -0.5, -0.5), 25))
        print("\nORDINARY")
        print(mp.nstr(ordinary(B_scaled, 0.5, 0.5), 25))
        print(mp.nstr(ordinary(B_scaled, 0.5, -0.5), 25))
        print(mp.nstr(ordinary(B_scaled, -0.5, 0.5), 25))
        print(mp.nstr(ordinary(B_scaled, -0.5, -0.5), 25))
        print("\nNEW")
        print(mp.nstr(New(B_scaled, 0.5, 0.5), 25))
        print(mp.nstr(New(B_scaled, 0.5, -0.5), 25))
        print(mp.nstr(New(B_scaled, -0.5, 0.5), 25))
        print(mp.nstr(New(B_scaled, -0.5, -0.5), 25))
        print("\nDIFFERENCE")
        print(mp.nstr(New(B_scaled, 0.5, 0.5) - ordinary(B_scaled, 0.5, 0.5), 40))
        print("\nCONTRIBUTION OF NEW PART")
        print((mp.mpf('0.5') * FINE_STRUCTURE_CONST**2 * B_scaled**3 * NUCLEAR_CHARGE**mp.mpf('-7/2') * 0.5))
        print((mp.mpf('0.5') * FINE_STRUCTURE_CONST ** 2 * B_scaled ** 3 * NUCLEAR_CHARGE ** mp.mpf('-7/2') * -0.5))


# GRAPHING
    Magnetic_Field = np.linspace(1, 50, 1024)
    Magnetic_Field_scaled = [mp.mpf(B) / mp.mpf('2.3505175e5') for B in Magnetic_Field]

    spin_combinations = [(0.5, 0.5), (0.5, -0.5), (-0.5, 0.5), (-0.5, -0.5)]

    Energy_Ordinary = [[ordinary(B, e_spin, n_spin) for B in Magnetic_Field_scaled] for e_spin, n_spin in spin_combinations]
    Energy_New = [[New(B, e_spin, n_spin) for B in Magnetic_Field_scaled] for e_spin, n_spin in spin_combinations]
    Ratio_Values = [[ratio(B, e_spin, n_spin) for B in Magnetic_Field_scaled] for e_spin, n_spin in spin_combinations]

    Graph_Energy(Magnetic_Field, Energy_New, [lbl for lbl in Labels])
    Graph_Ratio(Magnetic_Field, Ratio_Values, Labels)

if __name__ == '__main__':
    main()
