#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include "eigen.h"


const int Faces = 2;
const int Samples = 1;
const int Width = 92;
const int Height = 112;
const std::string Data_path = "/Users/olga_andreyeva/Dropbox/Developer/eigenface/faces/";
const int N = Faces*Samples; // Total number of images
const int M = Width*Height; // Resolution


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
        int val = image->array[0][i];
        if (val < 0)
        {
            val = 0;
        }
        image_file << val << " ";
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
        std::vector< std::vector<double> > facearray;
        for (int sample = 0; sample < Samples; ++sample)
        {
            std::stringstream filename;
            filename << Data_path << "s" << face + 1 << "/" << sample  + 1 << ".pgm";
            std::ifstream image(filename.str().c_str());

            if (image.is_open()) {
                facearray.push_back(read_pgm(image));
                image.close();
            } else {
                std::cout << "Image was not opened.";
            }
        }
        /* Find mean */
        std::vector<double> mean;
        for (int i = 0; i < M; ++i)
        {
            double sum = 0;
            for (int j = 0; j < Samples; ++j)
            {
                sum += facearray[j][i];
            }
            mean.push_back(sum/Samples);
        }
        array.push_back(mean);
    }
    return array;
}
/****/

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
Matrix multiply_columns(const Matrix& a, const Matrix& b) {
    Matrix result = Matrix(b.rows, a.columns);
    for (int c = 0; c < a.columns; ++c)
    {
        Matrix a_column = a.getColumn(c);
        result.setColumn(c, b*a_column);
    }
    return result;
}
/****/


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

