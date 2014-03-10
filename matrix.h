#include <iostream>
#include <vector>

class Matrix
{
public:
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
    Matrix getColumn(int number_of_column) const;
    Matrix getRow(int number_of_row) const;
    void setColumn(int number_of_column, Matrix vector);
    friend Matrix operator* (const Matrix& a, const Matrix& b);
};
