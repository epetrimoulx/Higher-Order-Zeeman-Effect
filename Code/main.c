/**
 * @author Evan Petrimoulx
 * @date December 8th 2024
 * @brief Program to solve the inhomogeneous perturbation equation by the method of Frobenius for an arbitrary state and calculate the higher order Zeeman effect
 *
 * @param argc The number of arguments
 * @param argv array containing the principle quantum number, the angular momentum, and the nuclear charge.
**/

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "perturbation_equation.h"

// Define constants here
const double EULER_MASCHERONI = 0.577215664901532860606512090082;
const double PI = 3.141592653589793238462643383279502884;

int main(int argc, char* argv[]) {
  int principle_quantum_num;
  int angular_momentum;
  int nuclear_charge;

  double psi_1_one_over_r[25] = {0}; // 1/r perturbation solution (First Order)  
  double psi_1_r_squared[25] = {0}; // r^2 perturbation solutions. (First Order)
  double psi_1_r_squared_2[25] = {0};
  double factorial[25] = {0};
  double integrals[50] = {0};
  double psi_norm[25] = {0};
  double harmonic_series[25] = {0};
  double inhomogeneous_terms[25] = {0};
  double recursion_coefficients[25] = {0};
  double hydrogenic_wavefunction[25] = {0};
  double energy[25] = {0}; // Each index is a higher order energy correction. Ex. energy[2] is the second order energy
  

  // Assign values to n, l, and Z based on program arguments
  switch(argc) {
    case 2:
      principle_quantum_num = atoi(argv[1]);
      angular_momentum = 0;
      nuclear_charge = 1;
      break;
    case 3:
      principle_quantum_num = atoi(argv[1]);
      angular_momentum = atoi(argv[2]);
      nuclear_charge = 1;
      break;
    case 4:
      principle_quantum_num = atoi(argv[1]);
      angular_momentum = atoi(argv[2]);
      nuclear_charge = atoi(argv[3]);
      break;
   default:
     principle_quantum_num = 1;
     angular_momentum = 0;
     nuclear_charge = 1;
     break;
  }

  printf("Starting Calculation of (H\u2070 - E\u2070) |\u03C8\u00B9\u27E9 = (1/r - E\u00B9) |\u03C8\u2070\u27E9\n");
  printf("N = %d\tL = %d\tZ = %d\n", principle_quantum_num, angular_momentum, nuclear_charge);
  
  // Compute and store factorials
  calc_factorials(25, factorial);
 
  // Compute and store harmonic_series
  calc_harmonic_series(25, harmonic_series);

  // Compute and store basic integrals
  calc_integrals(50, nuclear_charge, 2, factorial, integrals);
  
  energy[0] = -pow(nuclear_charge, 2) / 2.0;

  // For V = -1/r, needs to be changed for a different perturbation
  energy[1] = -nuclear_charge; 
    
  // Calculate the radial hydrogenic wavefunction coefficients
  calc_hydrogenic_wavefunction(nuclear_charge, angular_momentum, principle_quantum_num, factorial, hydrogenic_wavefunction);
  
  /*
   * The Hamiltonian Operator term causes the LHS of the power series solution to the perturbation equation to have a lowest power of r^(j-2). The summation starts at j = 0, so the lowest power in our solution is r^-2. The inhomogeneous terms array thus has inhomogeneous_terms[0] = the inhomogeneous term corresponding to r^{-2}. This piece of code will need to be modified if the solution requires powers less than -2.
   *
   * Starting with a 1/r perturbation, there are two inhomogeneous terms: 
       - r^-1
       - r^0
  */

  inhomogeneous_terms[1] = 1.0;        // -1/r term
  inhomogeneous_terms[2] = energy[1];  // r^0 term
  
  // Calculate recursion coefficients
  calc_recursion_coefficients(principle_quantum_num, angular_momentum, nuclear_charge, inhomogeneous_terms, integrals, recursion_coefficients);

  // Calculate a_0
  calc_overlap(integrals, hydrogenic_wavefunction, recursion_coefficients, psi_1_one_over_r);
  
  // Calculate psi_1 for the 1/r potential
  calc_psi_1(recursion_coefficients, psi_1_one_over_r);

  calc_norm(principle_quantum_num, psi_1_one_over_r, hydrogenic_wavefunction, integrals, psi_norm);

  for(int i = 0; i < 25; i++) {
    psi_1_one_over_r[i] = psi_norm[i];
  }

  
  printf("The recursion coefficiens for \u00B9/\u1D63 at L = 0 are:\n");
  for(int i = 0; i < 25; i++) {
    printf("%f\n", psi_1_one_over_r[i]);
  }

  /*
   * We will now compute the r^2 corrections to first order for l = 0
  */

  // First we have to reset all populated arrays that require different coefficients
  for(int i = 0; i < 25; i++) {
    inhomogeneous_terms[i] = 0.0;
    recursion_coefficients[i] = 0.0;
    psi_norm[i] = 0.0;
  }

  // Set up new first order energies
  energy[1] = 3.0 / pow(nuclear_charge, 2);

  // Set new inhomogeneous terms corresponding to r^2 case
  inhomogeneous_terms[2] = energy[1]; // r^0 term
  inhomogeneous_terms[4] = -1.0;      // r^2 term
  
  // Calculate the new recursion coefficients
  calc_recursion_coefficients(principle_quantum_num, angular_momentum, nuclear_charge, inhomogeneous_terms, integrals, recursion_coefficients);

  // Calculate a_0
  calc_overlap(integrals, hydrogenic_wavefunction, recursion_coefficients, psi_1_r_squared);
  
  // Calculate Psi_1
  calc_psi_1(recursion_coefficients, psi_1_r_squared);
  calc_norm(principle_quantum_num, psi_1_r_squared, hydrogenic_wavefunction, integrals, psi_norm);

  printf("The recursion coefficiens for r\u00B2 at L = 0 are:\n");
  for(int i = 0; i < 25; i++) {
    psi_1_r_squared[i] = psi_norm[i];
    printf("%f\n", psi_1_r_squared[i]);
  }

  /*
   * We will now compute the r^2 corrections to first order for l = 2
  */

  for(int i = 0; i < 25; i++) {
    inhomogeneous_terms[i] = 0.0;
    recursion_coefficients[i] = 0.0;
  }

  angular_momentum = 2;
  energy[1] = 0.0;

  inhomogeneous_terms[2] = energy[1];
  inhomogeneous_terms[4] = -1.0;  

  calc_recursion_coefficients(principle_quantum_num, angular_momentum, nuclear_charge, inhomogeneous_terms, integrals, recursion_coefficients);

  calc_psi_1(recursion_coefficients, psi_1_r_squared_2);

  calc_norm(principle_quantum_num, psi_1_r_squared_2, hydrogenic_wavefunction, integrals, psi_norm);

  printf("The recursion coefficiens for r\u00B2 at L = 2 are:\n");
  for(int i = 0; i < 25; i++) {
    psi_1_r_squared_2[i] = psi_norm[i];
    printf("%f\n", psi_1_r_squared_2[i]);
  }

  // NOT DONE BUT LOOKS PROMISING
  double one_over_r_matrix_element = calc_matrix_elements(psi_1_r_squared, hydrogenic_wavefunction, -1, -1, integrals);
  printf("\nRESULT r^2:\n%f\n", one_over_r_matrix_element);

  double r_squared_matrix_element = calc_matrix_elements(psi_1_one_over_r, hydrogenic_wavefunction, 1, 2, integrals);
  printf("\nRESULT 1/r:\n%f\n", r_squared_matrix_element);
}