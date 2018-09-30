//Author: Raphael Wagner 20.09.2018

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono> //used to measure runtime
#include "Jacobi_method.h"

void get_lowest_eigenvalue(double **A, double &lowest_eigenvalue, int N);
void fill_matrices(double **A, double **R, double potential(double),double h, int N);
double potential(double x);

using namespace std;
using namespace std::chrono;

int main()
{
    int gridpoints = 3;
    int gridpointsarray[] = {50,100,200};
    int max_values = 3;
    double rho_max_array[] = {2.0,4.0,8.0};
    double lowest_eigenvalue = 0;


    for (int n = 0; n < max_values; ++n) {
        double rho_max = rho_max_array[n];
        cout <<  "rho_max=" << rho_max << ":" << endl;

        for (int m = 0; m < gridpoints; ++m) {

            //set calculation parameter

            int N = gridpointsarray[m];
            double h = rho_max / double(N+1); //stepsize

            double tolerance = pow(10, -16);
            int rotations = 0;

            //initialize arrays

            auto **R = new double *[N]; //identity matrix to store Eigenvectors
            for (int i = 0; i < N; ++i) {
                R[i] = new double[N];
            }

            auto **A = new double *[N];
            for (int i = 0; i < N; ++i) {
                A[i] = new double[N];
            }

            //set matrix elements

            fill_matrices(A,R,potential,h,N);

            //Perform Jacobi Rotation

            //start time
            auto start = high_resolution_clock::now();

            Jacobi_EV(tolerance,rotations,A,R,N);

            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<nanoseconds>(stop - start);
            auto time = duration.count();

            get_lowest_eigenvalue(A,lowest_eigenvalue,N);

            //print results

            cout << "gridpoints:" << N << endl;
            cout << "lowest eigenvalue:" << lowest_eigenvalue << endl;
            cout << "runtime / seconds:" << time*pow(10,-9) << endl;
            cout << "rotations:" << rotations << endl;
            cout << endl;

            //destruct arrays

            for (int i = 0; i < N; ++i) {
                delete[] R[i];

            }
            delete[] R;

            for (int i = 0; i < N; ++i) {
                delete[] A[i];

            }
            delete[] A;
        }
    }
    return 0;
}

void get_lowest_eigenvalue(double **A, double &lowest_eigenvalue, int N)
{
    lowest_eigenvalue = A[0][0];

    for (int i = 0; i < N; ++i)
    {
        if (lowest_eigenvalue > A[i][i])
            lowest_eigenvalue = A[i][i];
    }
}

void fill_matrices(double **A, double **R, double potential(double), double h, int N)
{
    double hh = h * h;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; ++j)
        {
            if (i == j)
            {
                R[i][j] = 1;
                A[i][j] = 2.0/hh + potential(j*h+h);
            } else
            {
                R[i][j] = 0;
                A[i][j] = 0;
            }
            if ((i == (j + 1)) or (i == (j - 1)))
            {
                A[i][j] = -1.0/hh;
            }
        }
    }
}

double potential(double x)
{
    double y;
    y = x*x;
    return y;
}