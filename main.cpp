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
const std::string Data_path = "/Users/olga/Dropbox/Developer/eigenface/faces/";
const int N = Faces*Samples; // Total number of images
const int M = Width*Height; // Resolution


struct Matrix
{
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

};


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


Matrix transpose(Matrix *m) {
    Matrix result;
    result.rows = m->columns;
    result.columns = m->rows;
    for (int r = 0; r < result.rows; ++r)
    {
        std::vector<double> row;
        for (int c = 0; c < result.columns; ++c)
        {
            row.push_back(m->array[c][r]);
        }
        result.array.push_back(row);
    }
    return result;
}


Matrix mean_column(Matrix *m) {
    Matrix result;
    result.rows = m->rows;
    result.columns = 1;
    for (int r = 0; r < result.rows; ++r)
    {
        int sum = 0;
        std::vector<double> row;
        for (int c = 0; c < m->columns; ++c)
        {
            sum += m->array[r][c];
        }
        row.push_back(sum/m->columns);
        result.array.push_back(row);
    }
    return result;
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


Matrix multiply(Matrix *a, Matrix *b) {
    Matrix result;
    if (a->columns != b->rows)
    {
        std::cout << "Dimensions do not match" << std::endl;
        return result;
    }
    result.rows = a->rows;
    result.columns = b->columns;
    for (int r = 0; r < result.rows; r++)
    {
        std::vector<double> row (result.columns, 0);
        for (int c = 0; c < result.columns; ++c)
        {
            for (int k = 0; k < a->columns; ++k)
            {
                row[c] += a->array[r][k]*b->array[k][c];
            }
        }
        result.array.push_back(row);
    }
    return result;
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


Matrix multiply_columns(Matrix* a, Matrix*b) {
    Matrix result;
    Matrix at = transpose(a);
    result.rows = b->columns;
    result.columns = b->rows;
    for (int r = 0; r < result.rows; ++r)
    {
        Matrix eigenvector;
        eigenvector.rows = 1;
        eigenvector.columns = at.columns;
        eigenvector.array.push_back(at.array[r]);
        Matrix eigenvector_t = transpose(&eigenvector);
        eigenvector = multiply(b, &eigenvector_t);
        eigenvector = transpose(&eigenvector);
        result.array.push_back(eigenvector.array[0]);
    }
    return transpose(&result);
}


int main(int argc, const char * argv[])
{
    // First create matrix with images as rows
    Matrix images = Matrix(N, M, read_training_data());
    // Then transpose (images as columns)
    Matrix images_as_columns = transpose(&images);
    // Normalize images by subtracting the mean
    Matrix mean_image = mean_column(&images_as_columns);
    subtract_from_columns(&images_as_columns, &mean_image);

    // Transposed covariance matrix
    Matrix cov = multiply(&images, &images_as_columns);
    Matrix eigenvectors_of_cov = eigensystem(&cov).second;
    Matrix eigenfaces = multiply_columns(&eigenvectors_of_cov, &images_as_columns);

    return 0;
}

