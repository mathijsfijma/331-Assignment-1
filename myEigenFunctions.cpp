/*======================================================================================================

 This file is where you'll put the source information for the functions you write.  Make sure it's 
 included in your project (shift-alt-A) and that the functions you add are declared in the header
 file, myEigenFunctions.h.  You can add any extra comments or checks that could be relevant to these
 functions as you need to.

======================================================================================================*/

#include <string>

#include "myEigenFunctions.h"

using namespace std;


double DotProduct(double** A, double** B, int n, int m)
{
  //
  //	This is a function that takes two matrices A and B of identical dimensions (n*m) and 
  //  calculates and returns their dot product.
  //
  double dot = 0.0;

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      dot += A[i][j] * B[i][j];
    }
  }
  return dot;
}

double DotProduct(double* A, double* B, int n)
{
  //
  //	This is a function that takes two vectors A and B of identical length (n) and 
  //  calculates and returns their dot product.
  //
  double dot = 0.0;

  for (int i = 0; i < n; i++)
  {
    dot += A[i] * B[i];
  }
  return dot;
}

double* DotProduct(double** A, double* v, int n)
{
  //
  //  This is a function that takes a nxn-matrix A and an n-dimensional vector v stores
  //  the product A.v at the original location of v
  //
  double* result = new double[n]; // pointer to result vector

  for (int i = 0; i < n; i++)
  {
    result[i] = 0.0; // initialize ith element of result v
    for (int j = 0; j < n; j++)
    {
      result[i] += A[i][j] * v[j];
    }
  }

  return result;
}

double** ReadData(string inputFileName, int n, int m)
{
  //
  //  This is a function that returns a pointer to a 56x286 matrix, which contains the reflectance spectra.
  //  The first 29 rows contain spectra from Arabica samples, the following 27 rows contain spectra from Robusta samples.
  //  The code for parsing the CSV file has been adapted from Bernardo Trindade
  //  https://waterprogramming.wordpress.com/2017/08/20/reading-csv-files-in-c
  //

  double** spectra = new double*[n];
  for (int sample = 0; sample < n; sample++)
  {
    spectra[sample] = new double[m];
  }

  vector<vector<double> > data;
  cout << inputFileName << endl;
  ifstream inputFile(inputFileName);
  int l = 0;
  while (inputFile)
  {
    string s;
    if (!getline(inputFile, s)) break;
    // cout << s << endl;
    if (l != 0) // ignore first line, which contains wavelengths
    {
      istringstream ss(s);
      vector<double> record;

      int column = 0;
      while (ss)
      {
        // cout << "Row " << l << " " << ss << endl;
        string line;
        if (!getline(ss, line, ',')) break;
        try
        {
          // cout << "Row " << l << " Column " << column << " line " << line << endl; 
          spectra[l - 1][column] = stod(line);
          column++;
        }
        catch (const std::invalid_argument e)
        {
          cout << "NaN found in file " << inputFileName << " line " << l << endl;
          e.what();
        }
      }
    }
    l++;
  }

  return spectra;
}

void print_matrix(double** A, int n, int m)
{
  //
  // print_matrix prints the n-by-m matrix A.
  //
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      std::cout << A[i][j] << '\t';
    }
    std::cout << '\n';
  }
}

Eigenpair power_method(double** A, double* v, int n, double tol)
{
  //
  // This function computes the largest eigenvalue and the corresponding eigenvector
  //

  Eigenpair eigenpair(n);

  // Code Here //
    // initializing stuff
  eigenpair.length = n;
  eigenpair.vector = v;
  bool is_converged = false;
  double total = 0;
  eigenpair.normalize();
  double num = INFINITY;

  // iterative power method in while loop:
  while (!is_converged) {

      double* t = DotProduct(A, eigenpair.vector, n);
      double max = eigenpair.vector[0];

      // finding index of largest element in t vector
      int index = 0;
      for (int i = 1; i < n; ++i) {

          // Change < to > if you want to find the smallest element
          if (max < eigenpair.vector[i]) {
              max = eigenpair.vector[i];
              index = i;
          }
      }
      // finding estimate for eigenvalue and normalising eigenvector accordingly
      eigenpair.vector = t;
      eigenpair.normalize();
      
      if ((t[index] / eigenpair.vector[index]) < 0) {
          eigenpair.value = -eigenpair.value;
      }

      // convergance checks
      if (eigenpair.value < num) {
          if ((eigenpair.value + tol) >= num) {
              is_converged = true;
             // eigenpair.value = eigenpair.value;
              
          }
          num = eigenpair.value;
      }
      else if (eigenpair.value > num) {
          if ((eigenpair.value - num) <= tol) {
              is_converged = true;
              //eigenpair.value = eigenpair.value;
          }
          num = eigenpair.value;
      }

  }

  return eigenpair;
}

void deflate(double** A, Eigenpair eigenpair)
{
  //
  // This procedure removes eigenpair.vector from transformation matrix A in place
  //

  // Code Here //
    double temp_eigenvalue = eigenpair.value;
    eigenpair.normalize();


    for (int i = 0; i < eigenpair.length; i++) {
        for (int j = 0; j < eigenpair.length; j++) {
            A[i][j] = A[i][j] - (temp_eigenvalue * eigenpair.vector[i] * eigenpair.vector[j]);
        }
    }

}



double** CenterMatrix(double** A, int n, int m)
{
  //
  //  CenterMatrix takes a n-by-m matrix A and subtracts the emperical mean of each column.
  //  The function returns the pointer to the centered matrix.
  //

  double** result = new double*[n];

  for (int row = 0; row < n; row++)
  {
    result[row] = new double[m];
  }

    // Code Here //
// finding emperical means of each column of A
  double* means = new double[n];
  for (int i = 0; i < n; i++) {
      double total = 0;
      for (int j = 0; j < m; j++) {
          total += A[i][j];
      }
      means[i] = total / n;
  }
  // subtracting emperical mean from each column
  for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
          A[j][i] = A[j][i] - means[i];
      }
  }
  result = A;

  return result;
}


double** CovarianceMatrix(double** A, int n, int m)
{
  //
  //  CovarianceMatrix takes a n-by-m matrix A and computes its m-by-m covariance matrix.
  //  The function returns the pointer to the covariance matrix
  //

  double** cov = new double*[m];
  double** cA = CenterMatrix(A, n, m);

  for (int i = 0; i < m; i++)
  {
    cov[i] = new double[m];
  }

    // Code Here //
  for (int i = 0; i < m; i++) {
      for (int j = 0; j < m; j++) {
          double total = 0;
          for (int k = 0; k < n; k++) {
              total += cA[k][i] * cA[k][j];
          }
          cov[i][j] = total / n;
          // cov[i][j] = cA[n][i] * cA[n][j] / n;
      }
  }

  return cov;
}
