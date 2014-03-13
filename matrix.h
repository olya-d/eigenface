#include <iostream>
#include <vector>


/*
A basic matrix class.
With this class it is possible to:
- create 2D matrices
- print matrices
- tranpose matrices
- fetch specific row as a new matrix
- multiply matrices
*/
class Matrix
{
public:
    /* Attributes */
    int rows;
    int columns;
    std::vector< std::vector<double> > array;
    /* Constructors */
    Matrix();
    Matrix(int number_of_rows, int number_of_columns);
    Matrix(int number_of_rows, int number_of_columns, std::vector< std::vector<double> > elements);
    /* Methods */
    void print();
    Matrix transpose();
    Matrix getRow(int number_of_row) const;
    friend Matrix operator* (const Matrix& a, const Matrix& b);
};
