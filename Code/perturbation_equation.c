# include <stdio.h>
# include <math.h>
# include <stdlib.h>

/**
  * @author Evan Petrimoulx
  * @date January 24 2025
  * @brief Calculates any arbitrary E^(n)
  *
  * @param hydrogenic_wavefunction <psi_0|
  * @param perturbation The perturbing term
  * @param psi_n |psi_n>
  * @param integrals the array of integrals
  * @param power_of_r The power of r in the perturbing term. Ensures the power of r in the integral is correct
  * @return E^(n)
**/
double calc_energy(double* hydrogenic_wavefunction, double perturbation, double* psi_n, double* integrals, int power_of_r) {
  double energy = 0.0;
  
  for(int i = 0; i < 25; i++) {
    for(int j = 0; j < 25; j++) {
      energy = energy + hydrogenic_wavefunction[i] * perturbation * psi_n[j] * integrals[i + j + power_of_r + 2];
    }
  }

  return energy;
}


/**
  * @author Evan Petrimoulx
  * @date December 8th 2024 
  * @brief simple factorial function that stores all factorials terms in an array
  *
  * @param n the highest number in the factorial calculation
  * @param factorial the factorial array
  * @return factorial array
*/
void calc_factorials(int n, double* factorial){

  // if n is less than 0, the function is undefined, set n = 1
  if (n <= 0) {
    printf("Cannot have negative factorials, setting n = 1\n");
    n = 1;
  }
  
  factorial[0] = 1;
   
  for(int i = 1; i < n; i++) {
    factorial[i] = factorial[i - 1] * i;
  }
}


/**
  * @author Evan Petrimoulx
  * @date December 8th 2024
  * @brief Simple function to calculate the harmonic series up to n terms
  *
  * @param n the highest term in the summation
  * @param harmonic_series the array which stores the harmonic_series terms
  * @return harmonic_series array 
*/
void calc_harmonic_series(int n, double* harmonic_series){
  harmonic_series[0] = 1.0;

  for(int i = 1; i < n; i++) {
    harmonic_series[i] = harmonic_series[i - 1] + 1.0 / (i - 1);
  }
}


/**
  * @author Evan Petrimoulx
  * @date December 8th 2024 
  * @brief Calculates and stores integrals of the form  ∫ rʲ e⁻ᵅᶻʳ dr
  *
  * @param n the number of integrals to calculate. is also j in the mathematical expression.
  * @param nuclear_charge the nuclear charge (Z)
  * @param alpha the other varying parameter in the integral.
  * @param factorial the factorial array containing pre-computed factorials
  * @param integrals the array of integrals.
  * @return array of computed integrals
*/
void calc_integrals(int n, int nuclear_charge, int alpha, double* factorial, double* integrals) {
  
  integrals[0] = 1.0 / (alpha * nuclear_charge);
  
  for(int i = 1; i < n; i++) {
    integrals[i] =  factorial[i] / pow((alpha * nuclear_charge), i + 1);
  }
}

/**
  * @author Evan Petrimoulx
  * @date January 01 2025
  * @brief function to calculate the hydrogenic wavefunction coefficients for Psi_0, ignoring multiplying factors of r
  *
  * @param angular_momentum The angular momentum of the wavefunction
  * @param principle_quantum_num the principle quantum number n 
  * @param factorial an array of stored factorials
  * @param hydrogenic_wavefunction the array of hydrogenic_wavefunctions
  *
  * @note Calculates the Radial Wavefunction R_nl(r) = 2Z/n^2 sqrt(Z(n-l-1)! / (n+l)!) (2Zr/n)^l exp(-Zr/n) L^(2l+1)_(n-l-1) (2Zr/n). Where L is the Laguerre Polynomial. This function calculates everything in this equation but leaves out all multiplying factors of r.
*/
void calc_hydrogenic_wavefunction(int nuclear_charge, int angular_momentum, int principle_quantum_num, double* factorial, double* hydrogenic_wavefunction){
  double term;
  
  term = 2.0 * (double) nuclear_charge / pow(principle_quantum_num, 2) * pow(2.0 * (double) nuclear_charge / (double) principle_quantum_num, angular_momentum) * sqrt((double)nuclear_charge * factorial[principle_quantum_num - angular_momentum - 1] / factorial[principle_quantum_num + angular_momentum]);

  hydrogenic_wavefunction[0] = term;

  for(int i = 1; i < principle_quantum_num - angular_momentum - 1; i++){
    hydrogenic_wavefunction[i] = term * factorial[principle_quantum_num + angular_momentum + 1] / ((factorial[i] * factorial[principle_quantum_num - angular_momentum - 1 - i] * factorial[2 * angular_momentum + 2 + i]));   
  }
}

/**
  * @author Evan Petrimoulx
  * @date January 12th 2024
  * @brief function to calculate the norm of the hydrogenic wavefunction.
  *
  * @param principle_quantum_num the principle quantum number (n)
  * @param angular_momentum the angular momentum quantum number (l)
  * @param hydrogenic_wavefunction an array containing all of the coefficients of the hydrogenic radial wavefunction R_nl(r)
  * @param integrals an array containing all computed integrals needed for calculation
  * @param norm the norm that is returned after the calculation is completed.
  * @return the norm of the hydrogenic_wavefunction
*/
void calc_norm(int principle_quantum_num, double* psi_n, double* hydrogenic_wavefunction, double* integrals, double* psi_norm) {
  double inner_product = 0;
  double hydrogenic_norm = 0;
  double norm = 0;

  // Calculate <\psi_n | \psi_0> and <\psi_0|\psi_0>
  for(int i = 0; i < 25; i++){
    for(int j = 0; j < 25; j++) {
      inner_product = inner_product + psi_n[i] * hydrogenic_wavefunction[j] * integrals[i + j + 2];
      hydrogenic_norm = hydrogenic_norm + hydrogenic_wavefunction[i] * hydrogenic_wavefunction[j] * integrals[i + j + 2];
    }
  }

  double projection_coeff = (hydrogenic_norm != 0) ? (inner_product / hydrogenic_norm) : 0.0;

  for(int i = 0; i < 25; i++) {
    psi_norm[i] = psi_n[i] - projection_coeff * hydrogenic_wavefunction[i];
  }

  // **Compute the norm of psi_norm**
  for(int i = 0; i < 25; i++) {
    norm += psi_norm[i] * psi_norm[i];
  }

  norm = sqrt(norm);

  // **Explicitly normalize psi_norm**
  if (norm > 1e-10) {  // Avoid division by zero
    for(int i = 0; i < 25; i++) {
      psi_norm[i] /= norm;
    }
  }

  printf("Norm: %f\n", norm);
}

/**
  * @author Evan Petrimoulx 
  * @date January 12th 2024
  * @brief Function to calculate the recursion coefficients a_j for the power series solution to the perturbation equation
  *
  * @param principle_quantum_num the principle quantum number (n)
  * @param angular_momentum the angular momentum quantum number (l)
  * @param nuclear_charge The nuclear charge of the system. This is the quantum number (Z)
  * @param hydrogenic_wavefunction the array containing the coefficients of the hydrogenic radial wavefunction R_nl(r)
  * @param inhomogeneous_terms the array containing the inhomogeneous terms in the perturbation equation
  * @param integrals the array of integrals used to calculate the orthogonality condition
  * @param recursion_coefficients The resultant array with the calculated a_j values
  * @return the recursion coefficients
**/
void calc_recursion_coefficients(int principle_quantum_num, int angular_momentum, int nuclear_charge, double* inhomogeneous_terms, double* integrals, double* recursion_coefficients){
  
  // Calculate the terms in the series, look for terminating case
  recursion_coefficients[0] = inhomogeneous_terms[0] / (0.5 * (angular_momentum * (angular_momentum + 1)));
  recursion_coefficients[angular_momentum] = 0.0;
  
  for(int i = 1; i < 25; i++) {
    if(i == angular_momentum) {
      recursion_coefficients[i + 2] = 0.0;
      recursion_coefficients[i + 1] = inhomogeneous_terms[i + 2] / (nuclear_charge * (i + 1));
      recursion_coefficients[i] = (inhomogeneous_terms[i + 1] + 0.5 * ((i + 1) * (i + 2) - angular_momentum * (angular_momentum + 1)) * recursion_coefficients[i + 1]) / (nuclear_charge * i);

      i+=2;
      continue;
    }
    recursion_coefficients[i] = (-nuclear_charge * (i - 1) * recursion_coefficients[i - 1] + inhomogeneous_terms[i]) / (-0.5 * (i * (i + 1) - angular_momentum * (angular_momentum + 1)));
  }
}


/**
  * @author Evan Petrimoulx
  * @date January 17th 2024
  * @brief Function to calculate psi_1 given the recursion coefficients and psi_0
  *
  * @param integrals The array of calculated integrals
  * @param hydrogenic_wavefunction the array containing the coefficients of the radial hydrogenic wavefunction
  * @param recursion_coefficients the array containing the calculated recurstion coefficients of the power series expansion for psi_1 from the frobenius method
  * @param psi_1 the resultant array containing the final psi_1
  * @return psi_1
**/
void calc_overlap(double * integrals, double* hydrogenic_wavefunction, double* recursion_coefficients, double* psi_1) {
  double overlap_integral = 0;
  
  // Calculate <psi_0 | psi_1>
  for(int i = 0; i < 25; i++) {
    for(int j = 0; j < 25; j++) {
      overlap_integral = overlap_integral + recursion_coefficients[i] * hydrogenic_wavefunction[j] * integrals[i + j + 2];
    }
  }

  // Calculate psi_1
  for(int i = 0; i < 25; i++) {
    psi_1[i] = -overlap_integral * hydrogenic_wavefunction[i];
  }
}

/**
  * @author Evan Petrimoulx
  * @date January 24 2025
  * @brief Combines all recursion_coefficients and gets the final result
  *
  * @param psi_1 The final result
  * @param recursion_coefficients Array of all calculated recursion coefficients
  * @return psi_1
**/
void calc_psi_1(double* recursion_coefficients, double* psi_1){
  // Combine a_0 with the other a_j terms to get the final answer
  for(int i = 0; i < 25; i++) {
    psi_1[i] = psi_1[i] + recursion_coefficients[i];
  }
}