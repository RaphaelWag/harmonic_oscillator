//Documentations of functions are in the header

#include <iostream>
#include <cmath>
#include <random> //used to randomly check orthogonality of eigenvectors
#include "Jacobi_method.h"

using namespace std;

void Jacobi_rotate(double **A, double **R, double &maxvalue, int l, int k, int N)
{
    double tau = 0;
    double tan = 0;
    double cos = 0;
    double sin = 0;

    double A_kk, A_ll, A_ik, A_il, R_ik, R_il;

    if (maxvalue != 0) {

        //get trig. functions for rotation

        tau = (A[l][l] - A[k][k]) / (2 * A[k][l]);

        if (tau >= 0) {
            tan = 1.0 / (tau + sqrt(1.0 + tau * tau));
        } else {
            tan = -1.0 / (-tau + sqrt(1.0 + tau * tau));
        }

        cos = 1.0 / sqrt(1 + tan * tan);
        sin = cos * tan;

        if (tan == 0) {
            cos = 1.0;
            sin = 0.0;
        }

        A_kk = A[k][k];
        A_ll = A[l][l];
        A[k][k] = cos * cos * A_kk - 2.0 * cos * sin * A[k][l] + sin * sin * A_ll;
        A[l][l] = sin * sin * A_kk + 2.0 * cos * sin * A[k][l] + cos * cos * A_ll;
        A[k][l] = 0.0;
        A[l][k] = 0.0;

        for (int i = 0; i < N; i++) {
            if ((i != k) and (i != l)) {
                A_ik = A[i][k];
                A_il = A[i][l];
                A[i][k] = cos * A_ik - sin * A_il;
                A[k][i] = A[i][k];
                A[i][l] = cos * A_il + sin * A_ik;
                A[l][i] = A[i][l];
            }
            // And finally the new eigenvectors
            R_ik = R[i][k];
            R_il = R[i][l];
            R[i][k] = cos * R_ik - sin * R_il;

            R[i][l] = cos * R_il + sin * R_ik;
        }
    } else {
        cout << "search for biggest matrix element in Jacobi rotation has failed" << endl;
    }

}
void Jacobi_test(double **R, double tolerance, int N)
{
    int a = random_uniform(N);
    int b = random_uniform(N); // select 2 random eigenvectors and test orthonormality

    double S = 0;

    for (int i = 0; i < N; ++i)
    {
        S += R[i][a]*R[i][b];
    }
    if((a==b)and(abs(S-1)>tolerance)) {
        cout << "Jacobi_rotation_unittest() failed" << endl;
    }
    if((a!=b)and(abs(S)>tolerance)) {
        cout << "Jacobi_rotation_unittest() failed" << endl;;
    }
}



void get_max_element(double **A, double &maxvalue, int &l, int &k,int N)
{
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            if (abs(A[i][j]) > maxvalue) {
                maxvalue = abs(A[i][j]);
                l = i;
                k = j;
            }
        }
    }
}
int random_uniform(int N)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0,N-1); //N-1 because max indices of array is N-1
    return dis(gen);
}

void get_max_element_unittest(int N)
{
    double maxvalue =-1;
    int K;
    int L;

    //set up test matrix
    auto **A = new double *[N];
    for (int i = 0; i < N; ++i) {
        A[i] = new double[N];
    }

    //set all elements to 0

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i][j]=0;

        }
    }

    //set one random element nonzero
    int l = random_uniform(N);
    int k = random_uniform(N);

    if(k==l) k++; //we dont want a diagonal element to be the biggest
    if(k==N) k=k-2; //k should not be bigger than max. matrix indices

    A[l][k] = 1;

    get_max_element(A,maxvalue,L,K,N);

    if((k==K)and(l==L)) {
        cout << "get_max_element_unittest passed" << endl;
    }
    else{
        cout << "get_max_element_unittest failed" << endl;
    }

    //delete arrays

    for (int i = 0; i < N; ++i) {
        delete[] A[i];
    }
    delete[] A;
}

void Jacobi_EV(double tolerance, int &rotations, double **A, double **R, int N) {
    double maxvalue = tolerance + tolerance / 10.0;
    int l, k;

    while (maxvalue > tolerance) {
        maxvalue = 0;
        get_max_element(A, maxvalue, l, k, N);
        Jacobi_rotate(A, R, maxvalue, l, k, N);
        rotations++;
    }
}

void Jacobi_EV_unittest(double tolerance, double tolerance_2, int &rotations, double **A, double **R, int N)
{
    double maxvalue = tolerance + tolerance/10.0;
    int l,k;

    while (maxvalue > tolerance) {
        maxvalue = 0;
        get_max_element(A,maxvalue,l,k,N);
        Jacobi_rotate(A, R,maxvalue,l,k,N);
        rotations++;
        if (rotations % 100 == 0) { //check every 100 roations if eigenvectros are orthonormal
            Jacobi_test(R, tolerance_2, N);
        }
    }
    Jacobi_test(R,tolerance_2,N); // doing a final test
}
