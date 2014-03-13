#include <iostream>
#include <vector>
#include "matrix.h"

/*
Method to find the eigenvalues and eigenvectors of a symmetric matrix.
Returns the pair: first - vector of eigenvalues, second: matrix that has eigenvectors as columns.
Eigenvalues are sorted in the descending order, eigenvectors are sorted according to the eigenvalues in the descending order.
*/
std::pair<std::vector<double>, Matrix> eigensystem(Matrix *m);