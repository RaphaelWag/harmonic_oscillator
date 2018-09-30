#ifndef PROJECT_2_CPP_JACOBI_METHOD_H
#define PROJECT_2_CPP_JACOBI_METHOD_H

#endif //PROJECT_2_CPP_JACOBI_METHOD_H

/*
 * This Library provides all functions to solve an eigenvalue problem
 * including unit tests for the working functions
 *
 * Jacobi_rotate(double **A, double **R, double &maxvalue, int l, int k, int N)
 * This function performs a single rotation for Jacobis method on the matrix A.
 * The Matrix R should be initialized as an Identity matrix. The rotation gets
 * also applied to R so at the end this matrix contains the eigenvectors
 *
 * get_max_element(double **A, double &maxvalue, int &l, int &k,int N)
 * This function searches for the maximum off diagonal element of
 * a symmetric matrix and stores the information about the element
 * position in l and k.
 *
 * get_max_element_unittest(double maxvalue)
 * Unit test for function get_max_element
 *
 * Jacobi_rotation_unittest(double **R, double tolerance, int N)
 * This unit test randomly checks if the eingevectors stay orthonormal after
 * the rotation.
 *
 * random_uniform(int N)
 * This function generates random numbers from an uniform distribution in the range
 * from 0 to N-1
 *
 * Jacobi_EV(double tolerance,int &rotations, double **A, double **R, int N)
 * Performs Jacobi rotations until all off-diagonal elements are below the tolerance
 *
 */

void Jacobi_rotate(double **A, double **R, double &maxvalue, int l, int k, int N);
void get_max_element(double **A, double &maxvalue, int &l, int &k,int N);
void get_max_element_unittest(int N);
void Jacobi_test(double **R, double tolerance, int N, bool &test);
int random_uniform(int N);
void Jacobi_EV(double tolerance,int &rotations, double **A, double **R, int N);
void Jacobi_EV_unittest(double tolerance, double tolerance_2, int &rotations, double **A, double **R, int N);