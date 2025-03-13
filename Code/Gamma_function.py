import scipy.special
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import math

from sympy.core.facts import FactKB


def main() -> None:
    z = np.linspace(-5 - 5j, 15 + 15j, 512)
    x = np.linspace(0, 15, 15)
    Gamma_fn = scipy.special.gamma(z)  # Full complex Gamma(z)

    Factorial_fn = np.zeros(15)
    for i in range(0, 15):
        Factorial_fn[i] = math.factorial(i)

    print(Factorial_fn)

    # Seaborn and Matplotlib settings
    sns.set_style("whitegrid")
    plt.figure(figsize=(12, 8))

    # Plot Real and Imaginary parts of Gamma(z)
    plt.plot(z.real, Gamma_fn.real, color='blue', linewidth=2, label=r"$\Gamma(x)$")
    plt.plot(x, Factorial_fn, color = 'red', linewidth=2, label=r"$x!$")

    plt.yscale("log")
    # Titles and Labels
    plt.title("Gamma Function", fontsize=16, fontweight='bold')
    plt.xlabel("x", fontsize=14)
    plt.ylabel(r"$\Gamma(x)$", fontsize=14)

    # Legend
    plt.legend(fontsize=12)
    plt.grid(True, linestyle="--", linewidth=0.6)
    plt.show()
    plt.savefig("../Thesis/Gamma_function.pgf")

if __name__ == '__main__':
    main()
