#pragma once

void calc_psi_1(double* recursion_coefficients, double* psi_1);
void calc_factorials(int n, double* factorial);
void calc_harmonic_series(int n, double* harmonic_series);
void calc_integrals(int n, int nuclear_charge, int alpha, double* factorial, double* integrals);
void calc_overlap(double * integrals, double* hydrogenic_wavefunction, double* recursion_coefficients, double* psi_1);
void calc_hydrogenic_wavefunction(int nuclear_charge, int angular_momentum, int principle_quantum_num, double* factorial, double* hydrogenic_wavefunction);
void calc_recursion_coefficients(int angular_momentum, int nuclear_charge, double *inhomogeneous_terms, double *recursion_coefficients);
void calc_norm(double *psi_n, double *hydrogenic_wavefunction, double *integrals, double *psi_norm);
double calc_energy(const double* hydrogenic_wavefunction, double perturbation, const double* psi_n, const double* integrals, int power_of_r);
double calc_matrix_elements(double* psi_1, double* psi_2, double perturbation, int power_of_r, double* integrals);