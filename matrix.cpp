#include "matrix.h"
#include <iomanip>


Matrix::Matrix()
{
    rows = 0;
    columns = 0;
};
/* Constructors */
Matrix::Matrix(int number_of_rows, int number_of_columns)
{
    rows = number_of_rows;
    columns = number_of_columns;
    for (int r = 0; r < rows; ++r)
    {
        std::vector<double> row (columns, 0);
        array.push_back(row);
    }
}
Matrix::Matrix(int number_of_rows, int number_of_columns, std::vector< std::vector<double> > elements)
{
    rows = number_of_rows;
    columns = number_of_columns;
    if (elements.size() != number_of_rows)
    {
        std::cout << "Warning: number of rows does not match elements." << std::endl;
    }
    if (rows > 0 && elements[0].size() != number_of_columns)
    {
        std::cout << "Warning: number of columns does not match elements." << std::endl;
    }
    array = elements;
}
/* Methods */
void Matrix::print()
{
    std::cout << std::fixed;
    std::cout << std::setprecision(2);
    std::cout << "[" << std::endl;
    for (int r = 0; r < rows; ++r)
    {
        std::cout << "[";
        for (int c = 0; c < columns; ++c)
        {
            std::cout << array[r][c];
            if (c != columns - 1)
            {
                std::cout << ", ";
            }
        }
        std::cout << "]";
        if (r != rows - 1)
        {
            std::cout << "," << std::endl;
        }
    }
    std::cout << std::endl << "]" << std::endl;
}
Matrix Matrix::transpose()
{
    std::vector< std::vector<double> > result_array;
    for (int r = 0; r < columns; ++r)
    {
        std::vector<double> row;
        for (int c = 0; c < rows; ++c)
        {
            row.push_back(array[c][r]);
        }
        result_array.push_back(row);
    }
    return Matrix(columns, rows, result_array);
}
Matrix Matrix::getRow(int number_of_row) const
{
    if (number_of_row < 0 || number_of_row >= rows)
    {
        std::cout << "Error: row number is out of range" << std::endl;
        return Matrix();
    }
    std::vector< std::vector<double> > row;
    row.push_back(array[number_of_row]);
    return Matrix(1, columns, row);
}
Matrix operator* (const Matrix& a, const Matrix& b)
{
    if (a.columns != b.rows)
    {
        std::cout << "Error: dimensions do not match" << std::endl;
        return Matrix();
    }
    Matrix result = Matrix(a.rows, b.columns);
    for (int r = 0; r < a.rows; r++)
    {
        for (int c = 0; c < b.columns; ++c)
        {
            for (int k = 0; k < a.columns; ++k)
            {
                result.array[r][c] += a.array[r][k]*b.array[k][c];
            }
        }
    }
    return result;
}