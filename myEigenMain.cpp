/*======================================================================================================

    testEigenFunctions.cpp
	
	This is the main file for the ENGSCI331 Eigenvectors module.  
	It may demonstrate some new C++ syntax and functions.

    *** There are some examples of "bad" programming in here (bits missing etc.,) that you should
    probably fix. However, this file should compile straight away without any errors. ***

    You should use this file to get you started on the assignment.  
	You're welcome to change whatever you'd like to as you go, this is only a starting point.
   
======================================================================================================*/

#define _CRT_SECURE_NO_DEPRECATE

#include <iostream>
#include <fstream>
#include <string>

// This is the header file for your functions.  Usual programming practise would be to use a *.cpp 
// file that has the same name (ie: myEigenFunctions.cpp) and include it as normal in your project.
// Inside this file you'll see some ideas for functions that you could use during this project.  I 
// suggest you plan out your code first to see what kind of functions you'll use repeatedly and then
// write them.  
#include "myEigenFunctions.h"

using namespace std;

#define PI 3.14159265358979323846
#define TOL 1E-6

bool AE(double a, double b, double tol = TOL)
{
  return (abs(a - b) < tol); // Compare two values with respect to tolerance TOL
}


int main(void)
{
  /*==================================================================================================
  
      PART 1: Power method
  
  ==================================================================================================*/
  cout << "============================================================\n";
  cout << "PART 1: Power method\n\n";
  // Defining local array size variables
  int n = 2;

  // Defining and allocating memory for the 1D array v
  double* v = nullptr;
  v = new double [n];

  // Assigning values to the array v
  v[0] = 1.0;
  v[1] = 1.0;

  // Defining and allocating memory for the 2D array A - dimensions n-by-n:
  double** A = nullptr;
  A = new double*[n];
  for (int i = 0; i < n; i++)
  {
    A[i] = new double [n];
  }

  // Assigning values to the array A
  A[0][0] = 4.0;
  A[0][1] = 6.0;
  A[1][0] = 1.0;
  A[1][1] = 3.0;

  // Initializing struct comprising eigenvalue, eigenvector pair - see myEigenFunction.h
  Eigenpair eigenpair(n);

  // Calling power_method once
  eigenpair = power_method(A, v, n, TOL);

  if (AE(eigenpair.value, 6.0))
  {
    cout << "OK: (2x2) Power Method computing largest eigenvalue\n";
  }
  else
  {
    cout << "ERROR: (2x2) Power Method computing largest eigenvalue\n";
    cout << "eigenvalue: " << eigenpair.value << ", expected eigenvalue: 6\n";
  }

  if (AE(eigenpair.vector[0], 0.94868417, sqrt(TOL)) && AE(eigenpair.vector[1], 0.31622515, sqrt(TOL)))
  {
    cout << "OK: (2x2) Power Method computing eigenvector of largest eigenvalue\n";
  }
  else
  {
    cout << "ERROR: (2x2) Power Method computing eigenvector of largest eigenvalue\n";
    cout << "eigenvector         (" << eigenpair.vector[0] << "," << eigenpair.vector[1] << ")\n";
    cout << "expected eigenvector(0.94868417 , 0.31622515)\n";
  }

  delete[] v;
  cout << "============================================================\n";
  /*==================================================================================================
  
      PART 2: Power Method and Deflate
  
  ==================================================================================================*/
  cout << "PART 2: Power Method and Deflate\n\n";
  // Defining and allocating memory for the 2D array C - dimensions nC-by-nC:
  int nC = 3;
  double** C = nullptr;
  C = new double*[nC];
  for (int i = 0; i < nC; i++)
  {
    C[i] = new double [nC];
  }

  // Assigning values from the C.txt file
  ifstream infile;
  infile.open("C.txt");
  for (int i = 0; i < nC; i++)
  {
    for (int j = 0; j < nC; j++)
    {
      infile >> C[i][j];
    }
  }
  infile.close();

  // Checking if file read properly.
  if (AE(C[0][0], 1.0) && AE(C[nC - 1][nC - 1], 2.0))
  {
    cout << "OK: Reading of the file C.txt" << endl;
  }
  else
  {
    cout << "ERROR: Reading of the file C.txt, make sure it is in the correct folder" << endl;
  }
  cout << "Test matrix C:" << endl;
  print_matrix(C, nC, nC);

  // Redefining and allocating memory for the 1D array v
  v = new double [nC];

  // Assigning values to the array v
  v[0] = 1.0;
  v[1] = 1.0;
  v[2] = 1.0;

  // Initializing struct comprising eigenvalue, eigenvector pair - see myEigenFunction.h
  Eigenpair eigenpairC(nC);

  // Calling power_method once
  eigenpairC = power_method(C, v, nC, TOL);

  if (AE(eigenpairC.value, 3.0))
  {
    cout << "OK: (3x3) Power Method computing largest eigenvalue\n";
  }
  else
  {
    cout << "ERROR: (3x3) Power Method computing largest eigenvalue\n";
    cout << "eigenvalue: " << eigenpairC.value << ", expected eigenvalue: 3\n";
  }

  if (AE(eigenpairC.vector[0], 0.70710678, sqrt(TOL)) &&
    AE(eigenpairC.vector[1], 0.70710678, sqrt(TOL)) &&
    AE(eigenpairC.vector[2], 0, sqrt(TOL)))
  {
    cout << "OK: (3x3) Power Method computing eigenvector of largest eigenvalue\n";
  }
  else
  {
    cout << "ERROR: (3x3) Power Method computing eigenvector of largest eigenvalue\n";
    cout << "eigenvector         (" << eigenpairC.vector[0] << "," << eigenpairC.vector[1] << "," << eigenpairC.vector[
      2] << ")\n";
    cout << "expected eigenvector(0.70710678 , 0.70710678,  0)\n";
  }

  cout << "Calling deflate\n";
  deflate(C, eigenpairC);
  if (AE(C[0][0], -0.5, cbrt(TOL)) && AE(C[0][1], 0.5, cbrt(TOL)) && AE(C[0][2], 0.0, cbrt(TOL)) &&
    AE(C[1][0], 0.5, cbrt(TOL)) && AE(C[1][1], -0.5, cbrt(TOL)) && AE(C[1][2], 0.0, cbrt(TOL)) &&
    AE(C[2][0], 0.0, cbrt(TOL)) && AE(C[2][1], 0.0, cbrt(TOL)) && AE(C[2][2], 2.0, cbrt(TOL))
  )
  {
    cout << "OK: deflate method computes deflated matrix correctly\n";
  }
  else
  {
    cout << "ERROR: deflate method does not compute deflated matrix correctly\n";
    print_matrix(C, nC, nC);
  }

  // Are these what you expect? Perhaps you should work out what they should be.
  cout << "Calling power_method for the second time\n";
  eigenpairC = power_method(C, v, nC, TOL);
  eigenpairC.print();

  cout << "Calling deflate for the second time\n";
  deflate(C, eigenpairC);
  print_matrix(C, nC, nC);

  cout << "Calling power_method for the third time\n";
  eigenpairC = power_method(C, v, nC, TOL);
  eigenpairC.print();
  cout << "============================================================\n";
  /*==================================================================================================
  
      PART 3: Center and Covariance
  
  ==================================================================================================*/
  cout << "PART 3: Center and Covariance\n\n";
  C = new double*[nC];
  for (int i = 0; i < nC; i++)
  {
    C[i] = new double [nC];
  }

  // Assigning values from the C.txt file
  infile.open("C.txt");
  for (int i = 0; i < nC; i++)
  {
    for (int j = 0; j < nC; j++)
    {
      infile >> C[i][j];
    }
  }
  infile.close();

  // Checking if file read properly.
  if (AE(C[0][0], 1.0) && AE(C[nC - 1][nC - 1], 2.0))
  {
    cout << "OK: Reading of the file C.txt" << endl;
  }
  else
  {
    cout << "ERROR: Reading of the file C.txt, make sure it is in the correct folder" << endl;
  }
  cout << "Test matrix C:" << endl;
  print_matrix(C, nC, nC);

  // Are these what you expect? Perhaps you should work out what they should be.
  cout << "Centering Matrix" << endl;
  double** CC = CenterMatrix(C, nC, nC);
  print_matrix(CC, nC, nC);

  cout << "Covariance Matrix" << endl;
  double** cov = CovarianceMatrix(C, nC, nC);
  print_matrix(cov, nC, nC);

  // Writing Cov Matrix to file
  ofstream cov_out_file;
  cov_out_file.open("cov.txt");
  for (int i = 0; i < nC; i++)
  {
    for (int j = 0; j < nC; j++)
    {
      cov_out_file << cov[i][j];
      cov_out_file << (j != nC - 1 ? "," : "\n");
    }
  }
  cov_out_file.close();

  cout << "============================================================\n";
  /*==================================================================================================
  
      PART 4: Spectra Data
  
  ==================================================================================================*/
  cout << "PART 4: Spectra Data\n\n";
  const int N = 56; // rows
  const int M = 286; // columns

  double** spectra = ReadData("DS19hH2_dk0_FTIR_Spectra_instant_coffee.csv", N, M);

  // Test reading of CSV file
  if (AE(spectra[0][0], 21.227620) && AE(spectra[N - 1][M - 1], 1.574679))
  {
    cout << "OK: Reading of spectra CSV file" << endl;
  }
  else
  {
    cout << "ERROR: Reading of spectra CSV file" << endl;
  }

  cov = CovarianceMatrix(spectra, N, M);

  // Code Here //
  double vector_guess[N];
  double eigenvalues[N];
  double* eigenvectors[N];
  
  //computing principle components, iteratively calling power_method and deflate
  for (int i = 0; i < N; i++) {
      
      Eigenpair Eigens = power_method(cov, vector_guess, (N), TOL);
      eigenvalues[i] = Eigens.value;
      eigenvectors[i] = Eigens.vector;
      deflate(cov, Eigens);
  }

  //writing out eigenvectors and eigenvalues (representing principle components) to csv
  std::ofstream PCA_FILE_EIGENVECTORS;
  PCA_FILE_EIGENVECTORS.open("PCA_EIGENVECTOR.csv");
  PCA_FILE_EIGENVECTORS << "PRINCIPAL COMPONENTS OF SPECTRA COVARIANCE MATRIX.\n";
  for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
          PCA_FILE_EIGENVECTORS << eigenvectors[i][j] << ",";
      }
      PCA_FILE_EIGENVECTORS << ".\n";
  }

  std::ofstream PCA_FILE_EIGENVALUES;
  PCA_FILE_EIGENVALUES.open("PCA_EIGENVALUES.csv");
  PCA_FILE_EIGENVALUES << "PRINCIPAL COMPONENTS OF SPECTRA COVARIANCE MATRIX.\n";
  for (int i = 0; i < N; i++) {
      PCA_FILE_EIGENVALUES << eigenvalues[i] << ",";
  }

// to do: ORTHOGONALITY CHECKS

  cout << "============================================================\n";
}
