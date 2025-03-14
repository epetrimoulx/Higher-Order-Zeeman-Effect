import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def main() -> None:
    # Define x values for plotting
    x_vals = np.linspace(-1, 1, 512)

    # Compute the first five Legendre polynomials
    x = sp.symbols('x')
    legendre_polynomials = [sp.legendre(n, x) for n in range(5)]
    legendre_funcs = [sp.lambdify(x, poly, 'numpy') for poly in legendre_polynomials]

    # Evaluate polynomials
    legendre_values = [func(x_vals) for func in legendre_funcs]

    # Plot setup
    sns.set_style("whitegrid")
    plt.figure(figsize=(12, 8), dpi=200)

    # Plot the polynomials
    colors = ["blue", "red", "green", "black", "purple"]
    for i in range(5):
        sns.lineplot(x=x_vals, y=legendre_values[i], color=colors[i], linewidth=2, label=rf"$P_{i}(x)$")

    # Titles and Labels
    plt.title("First Few Legendre Polynomials", fontsize=16, fontweight='bold')
    plt.xlabel("x", fontsize=14)
    plt.ylabel(r"$P_l(x)$", fontsize=14)

    # Legend
    plt.legend(fontsize=12, fancybox=True)

    # Save and Show
    plt.savefig("../Thesis/Legendre_polynomials.pgf")
    plt.show()


if __name__ == '__main__':
    main()
