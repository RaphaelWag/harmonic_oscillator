//Author: Raphael Wagner 20.09.2018

#include <iostream>
#include <fstream>
#include <cmath>
#include <random> //used to randomly check orthogonality of eigenvectors

void Jacobi_rotate(double **A, double **R, double &maxvalue, int l, int k, int N);
void get_max_element(double **A, double &maxvalue, int &l, int &k,int N);
void Jacobi_rotation_unittest(double **R, double tolerance, int N);
int random_uniform(int N);
void print_results_txt(double **results, int N);

using namespace std;

int main()
{
    //set calculation parameter

    int N = 101;
    double rho_max = 4.75;

    double h = rho_max/double(N); //stepsize

    double hh = h*h;

    double tolerance = pow(10,-16);
    double tolerance_2 = pow(10,-8);
    double maxvalue = 0;
    maxvalue = tolerance + tolerance/10.0;
    int l;
    int k;

    auto **R = new double *[N]; //identity matrix to store Eigenvectors
    for (int i = 0; i < N; ++i) {
        R[i] = new double[N];
    }

    auto **A = new double *[N];
    for (int i = 0; i < N; ++i) {
        A[i] = new double[N];
    }


    //set array values
    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            if (i == j)
            {
                R[i][j] = 1;
                A[i][j] = 2.0/hh + (i)*(i)*hh;
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

    //Perform Jacobi Rotation

    //start time

    while (maxvalue > tolerance) {
        maxvalue = 0;
        get_max_element(A,maxvalue,l,k,N);
        Jacobi_rotate(A, R,maxvalue,l,k,N);
    }
    Jacobi_rotation_unittest(R,tolerance_2,N); // doing a final unit test just to be sure


    for (int p = 0; p < N; ++p) {
        cout << A[p][p] << endl;
    }

    //print results for python

    print_results_txt(A,N);

    return 0;
}

void Jacobi_rotate(double **A, double **R, double &maxvalue, int l, int k, int N)
{
    //search for highest valued off diagonal element
    double tau = 0;
    double t = 0;
    double c = 0;
    double s = 0;

    if (maxvalue != 0) {


        //get trig. functions for rotation

        tau = (A[l][l] - A[k][k]) / (2.0 * A[k][l]);

        if (tau >= 0) {
            t = 1.0 / (tau + sqrt(1.0 + tau * tau));
        } else {
            t = -1.0 / (-tau + sqrt(1.0 + tau * tau));
        }

        c = 1.0 / sqrt(1 + t * t);
        s = c * t;

        if (t == 0) {
            c = 1.0;
            s = 0.0;
        }

        double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
        a_kk = A[k][k];
        a_ll = A[l][l];
        A[k][k] = c * c * a_kk - 2.0 * c * s * A[k][l] + s * s * a_ll;
        A[l][l] = s * s * a_kk + 2.0 * c * s * A[k][l] + c * c * a_ll;
        A[k][l] = 0.0; // hard-coding non-diagonal elements by hand
        A[l][k] = 0.0; // same here
        for (int i = 0; i < N; i++) {
            if ((i != k) and (i != l)) {
                a_ik = A[i][k];
                a_il = A[i][l];
                A[i][k] = c * a_ik - s * a_il;
                A[k][i] = A[i][k];
                A[i][l] = c * a_il + s * a_ik;
                A[l][i] = A[i][l];
            }
            // And finally the new eigenvectors
            r_ik = R[i][k];
            r_il = R[i][l];
            R[i][k] = c * r_ik - s * r_il;

            R[i][l] = c * r_il + s * r_ik;
        }
    } else {
        cout << "search for biggest matrix element in Jacobi rotation has failed" << endl;
    }

}

void print_runtime_txt(double *runtimearray, int *rotationarray,int *gridpoints, int N)
{
    string path = "D:\\Uni\\Lectures\\Computational Physics\\Project 2 Codes\\project_2_python_tests\\"; //set your path to your python project here to print results in python
    string content("runtime_rotations_");
    string txt(".txt");
    string file;
    file = path+content+txt;

    ofstream myfile(file);
    if (!myfile.is_open())
        cout << "Unable to open file";
    else {
        for (int k = 0; k < N; k++) {
            myfile << gridpoints[k] << " " << log10(runtimearray[k]) << " " << rotationarray[k] << " ";
        }
        myfile.close();
    }
}

void Jacobi_rotation_unittest(double **R,double tolerance, int N)
{
    int a = random_uniform(N);
    int b = random_uniform(N);

    double S = 0;

    for (int i = 0; i < N; ++i)
    {
        S += R[i][a]*R[i][b];
    }
    if((a==b)and(abs(S-1)>tolerance))
    {
        cout << "Eigenvector orthogonality test failed" << endl;
    }
    if((a!=b)and(abs(S)>tolerance))
    {
        cout << "Eigenvector orthogonality test failed" << endl;
    }
}



void get_max_element(double **A, double &maxvalue, int &l, int &k,int N)
{
    for (int i = 1; i < N; ++i) {
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
    uniform_int_distribution<> dis(0,N-1);
    return dis(gen);
}

void print_results_txt(double **results, int N)
{
    string path = "D:\\Uni\\Lectures\\Computational Physics\\Project 2 Codes\\one_electron_python\\"; //set your path to your python project here to print results in python
    string content("eigenvalues");
    string txt(".txt");
    string file;
    file = path+content+txt;

    ofstream myfile(file);
    if (!myfile.is_open())
        cout << "Unable to open file";
    else {
        for (int k = 0; k < N; k++) {
            myfile << results[k][k] << " ";
        }
        myfile.close();
    }
}