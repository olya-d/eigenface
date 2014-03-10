#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>


const int Faces = 2;
const int Samples = 1;
const int Width = 92;
const int Height = 112;
const std::string Data_path = "/Users/olga_andreyeva/Dropbox/Developer/eigenface/faces/";
const int N = Faces*Samples; // Total number of images
const int M = Width*Height; // Resolution


class Matrix
{
public:
    int rows;
    int columns;
    std::vector< std::vector<double> > array;
    Matrix()
    {
        rows = 0;
        columns = 0;
    };
    /* Constructors */
    Matrix(int number_of_rows, int number_of_columns)
    {
        rows = number_of_rows;
        columns = number_of_columns;
        for (int r = 0; r < rows; ++r)
        {
            std::vector<double> row (columns, 0);
            array.push_back(row);
        }
    }
    Matrix(int number_of_rows, int number_of_columns, std::vector< std::vector<double> > elements)
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
    void print()
    {
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
    Matrix transpose()
    {
        // Matrix result = Matrix(columns, rows);
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
    Matrix getColumn(int number_of_column) const
    {
        if (number_of_column < 0 || number_of_column >= columns)
        {
            std::cout << "Error: column number is out of range" << std::endl;
            return Matrix();
        }
        std::vector< std::vector<double> > column_values;
        for (int r = 0; r < rows; ++r)
        {
            std::vector<double> row;
            row.push_back(array[r][number_of_column]);
            column_values.push_back(row);
        }
        return Matrix(rows, 1, column_values);
    }
    void setColumn(int number_of_column, Matrix vector)
    {
        if (number_of_column < 0 || number_of_column >= columns)
        {
            std::cout << "Error: column number is out of range" << std::endl;
            return;
        }
        if (rows != vector.rows)
        {
            std::cout << "Error: number of rows does not match" << std::endl;
            return;
        }
        for (int r = 0; r < rows; ++r)
        {
            array[r][number_of_column] = vector.array[r][0];
        }
    }
    friend Matrix operator* (const Matrix& a, const Matrix& b)
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
};

/* Work with data */
std::vector<double> read_pgm(std::ifstream& file, int size=M) {
    std::vector<double> values;
    std::string line;
    getline(file, line); // Skip P2 line
    getline(file, line); // Skip width line
    getline(file, line); // Skip height line
    getline(file, line); // Skip max value line

    int val;
    while(file >> val)
    {
        values.push_back(val);
    }
    return values;
}

void write_pgm(std::string file, Matrix *image) {
    std::stringstream filename;
    filename << file;
    std::ofstream image_file(filename.str().c_str());
    image_file << "P2" << std::endl << Width << std::endl << Height << std::endl << "255" << std::endl;
    for (int i = 0; i < M; ++i)
    {
        image_file << image->array[i][0] << " ";
    }
    image_file.close();
}

std::vector< std::vector<double> > read_training_data() {
    /*
     Returns pointer to the NxM array a, s.t
     a[i][j] is jth value of ith image.
     */
    std::vector< std::vector<double> > array;
    for (int face = 0; face < Faces; ++face)
    {
        for (int sample = 0; sample < Samples; ++sample)
        {
            std::stringstream filename;
            filename << Data_path << "s" << face + 1 << "/" << sample  + 1 << ".pgm";
            std::ifstream image(filename.str().c_str());

            if (image.is_open()) {
                array.push_back(read_pgm(image));
                image.close();
            } else {
                std::cout << "Image was not opened.";
            }
        }
    }
    return array;
}

/* Work with matrices */
Matrix mean_column(const Matrix& m) {
    std::vector< std::vector<double> > array;
    for (int r = 0; r < m.rows; ++r)
    {
        double sum = 0;
        std::vector<double> row;
        for (int c = 0; c < m.columns; ++c)
        {
            sum += m.array[r][c];
        }
        row.push_back(sum/m.columns);
        array.push_back(row);
    }
    return Matrix(m.rows, 1, array);
}


void subtract_from_columns(Matrix *m, Matrix *v) {
    /*
     Subtracts std::vector v from each column of m.
     */
    for (int r = 0; r < m->rows; ++r)
    {
        for (int c = 0; c < m->columns; ++c)
        {
            m->array[r][c] -= v->array[r][0];
            if (m->array[r][c] < 0)
            {
                m->array[r][c] = 0;
            }
        }
    }
}


Matrix identity_matrix(int dimension)
{
    Matrix I;
    I.rows = dimension;
    I.columns = dimension;
    for (int r = 0; r < dimension; ++r)
    {
        std::vector<double> row;
        for (int c = 0; c < dimension; ++c)
        {
            if (c == r)
            {
                row.push_back(1);
            }
            else
            {
                row.push_back(0);
            }
        }
        I.array.push_back(row);
    }
    return I;
}


std::pair<std::vector<double>, Matrix> eigensystem(Matrix *m) {
    /*
     Using Jacobi eigenvalue method, returns pair:
     vector of eigenvalues and matrix, which columns are eigenvectors.
     */
    Matrix A; // will become diagonal
    A.rows = m->rows;
    A.columns = m->columns;
    A.array = m->array;
    Matrix P = identity_matrix(A.rows); // accumulates rotations
    int max_rotations = 5*A.rows*A.rows;
    for (int i = 0; i < A.rows; ++i)
    {
        std::vector<double> row;
        for (int j = 0; j < A.columns; ++j)
        {
            row.push_back(m->array[i][j]);
        }
        A.array.push_back(row);
    }

    std::vector<double> eigenvalues;

    for (int it = 0; it < max_rotations; ++it)
    {
        double max = 0, k, l;
        // Find the largest off-diagonal element
        for (int r = 0; r < A.rows - 1; ++r)
        {
            for (int c = r + 1; c < A.rows; ++c)
            {
                if (abs(A.array[r][c]) >= max)
                {
                    max = abs(A.array[r][c]);
                    k = r;
                    l = c;
                }
            }
        }
        if (max < 1.0e-12)
        {
            for (int i = 0; i < A.rows; ++i)
            {
                eigenvalues.push_back(A.array[i][i]);
            }
            // Normalize P
            for (int c = 0; c < P.columns; ++c)
            {
                double lenght = 0;
                for (int r = 0; r < P.rows; ++r)
                {
                    lenght += P.array[r][c]*P.array[r][c];
                }
                for (int r = 0; r < P.rows; ++r)
                {
                    P.array[r][c] = P.array[r][c]/lenght;
                }
            }

            std::pair<std::vector<double>, Matrix> result(eigenvalues, P);
            return result;
        }
        // Perform rotation
        double diff = A.array[l][l] - A.array[k][k];

        double t;
        if (abs(A.array[k][l]) < abs(diff)*1.0e-36)
        {
            t = A.array[k][l]/diff;
        }
        else
        {
            double phi = diff/(2.0*A.array[k][l]);
            t = 1.0/(abs(phi) + sqrt(phi*phi + 1.0));
            if (phi < 0)
            {
                t = -t;
            }
        }
        double c = 1.0/sqrt(t*t + 1.0);
        double s = t*c;
        double tau = s/(1.0 + c);
        double temp = A.array[k][l];
        A.array[k][l] = 0;
        A.array[k][k] = A.array[k][k] - t*temp;
        A.array[l][l] = A.array[l][l] + t*temp;

        for (int i = 0; i < k; ++i)
        {
            temp = A.array[i][k];
            A.array[i][k] = temp - s*(A.array[i][l] + tau*temp);
            A.array[i][l] = A.array[i][l] + s*(temp - tau*A.array[i][l]);
        }
        for (int i = k + 1; i < l; ++i)
        {
            temp = A.array[k][i];
            A.array[k][i] = temp - s*(A.array[i][l] + tau*A.array[k][i]);
            A.array[i][l] = A.array[i][l] + s*(temp - tau*A.array[i][l]);
        }
        for (int i = l + 1; i < A.rows; ++i)
        {
            temp = A.array[k][i];
            A.array[k][i] = temp - s*(A.array[l][i] + tau*temp);
            A.array[l][i] = A.array[l][i] + s*(temp - tau*A.array[l][i]);
        }
        for (int i = 0; i < A.rows; ++i)
        {
            temp = P.array[i][k];
            P.array[i][k] = temp - s*(P.array[i][l] + tau*P.array[i][k]);
            P.array[i][l] = P.array[i][l] + s*(temp - tau*P.array[i][l]);
        }
    }
    std::cout << "Jacobi method didn't converge." << std::endl;
    for (int i = 0; i < A.rows; ++i)
    {
        eigenvalues.push_back(A.array[i][i]);
    }
    std::pair<std::vector<double>, Matrix> result(eigenvalues, P);
    return result;
}


Matrix multiply_columns(const Matrix& a, const Matrix& b) {
    Matrix result = Matrix(b.rows, a.columns);
    for (int c = 0; c < a.columns; ++c)
    {
        Matrix a_column = a.getColumn(c);
        result.setColumn(c, b*a_column);
    }
    return result;
}


int main(int argc, const char * argv[])
{
    // First create matrix with images as rows
    Matrix images = Matrix(N, M, read_training_data());
    // Then transpose (images as columns)
    Matrix images_as_columns = images.transpose();
    // Normalize images by subtracting the mean
    Matrix mean_image = mean_column(images_as_columns);
    subtract_from_columns(&images_as_columns, &mean_image);
    images = images_as_columns.transpose();

    // Transposed covariance matrix
    Matrix cov = images*images_as_columns;
    cov.print();
    Matrix eigenvectors_of_cov = eigensystem(&cov).second;
    Matrix eigenfaces = multiply_columns(eigenvectors_of_cov, images_as_columns);

    return 0;
}

