import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pyx.unit import x_inch


def main() -> None:
    # Define x values for plotting
    x_vals = np.linspace(-5, 15, 512)

    laguerre_polynomials = np.zeros((5, 512))
    laguerre_polynomials[0] = 1 + 0*x_vals
    laguerre_polynomials[1] = -x_vals + 1
    laguerre_polynomials[2] = 0.5 * (x_vals**2 - 4*x_vals + 2)
    laguerre_polynomials[3] = 1.0/6.0 * (-x_vals**3 + 9*x_vals**2 - 18*x_vals + 6)
    laguerre_polynomials[4] = 1.0/24.0 * (x_vals**3 - 16*x_vals**3 + 72*x_vals**2 - 96*x_vals +24)


    # Plot setup
    sns.set_style("whitegrid")
    plt.figure(figsize=(12, 8), dpi=200)

    # Plot the polynomials
    colors = ["blue", "red", "green", "black", "purple"]
    for i in range(5):
        sns.lineplot(x=x_vals, y=laguerre_polynomials[i], color=colors[i], linewidth=2, label=rf"$L_{i}(x)$")

    # Titles and Labels
    plt.title("First Few Laguerre Polynomials", fontsize=16, fontweight='bold')
    plt.xlabel("x", fontsize=14)
    plt.ylabel(r"$L_n(x)$", fontsize=14)
    plt.ylim(-20, 20)
    plt.xlim(-5, 10)
    # Legend
    plt.legend(fontsize=12, fancybox=True)

    # Axes and Grid
    plt.axhline(0, color="black", linewidth=1, linestyle="--")
    plt.axvline(0, color="black", linewidth=1, linestyle="--")
    plt.grid(True, linestyle="--", linewidth=0.6)
    # Save and Show
    plt.savefig("../Thesis/Laguerre_polynomials.pgf")
    plt.show()


if __name__ == '__main__':
    main()
