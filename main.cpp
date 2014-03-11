#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include "eigen.h"


const int Faces = 40;
const int Samples = 5;
const int Width = 92;
const int Height = 112;
const int Eigenfaces = 27;
const std::string Data_path = "/Users/olga/Dropbox/Developer/eigenface/faces/";
const int N = Faces; // Total number of images
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
Matrix scale(Matrix m, double min = 0, double max = 255) {
    double m_min = m.array[0][0];
    double m_max = m.array[0][0];
    for (int r = 0; r < m.rows; ++r) {
        for (int c = 0; c < m.columns; ++c) {
            if (m.array[r][c] < m_min) {
                m_min = m.array[r][c];
            }
            if (m.array[r][c] > m_max) {
                m_max = m.array[r][c];
            }
        }
    }
    double old_range = m_max - m_min;
    double new_range = max - min;

    Matrix result;
    result.columns = m.columns;
    result.rows = m.rows;
    for (int r = 0; r < m.rows; ++r) {
        std::vector<double> row;
        for (int c = 0; c < m.columns; ++c) {
            row.push_back((m.array[r][c] - m_min)*new_range/old_range + min);
        }
        result.array.push_back(row);
    }
    return result;
}
/****/


int main(int argc, const char * argv[])
{
    /*
    A is a NxM matrix, s.t. A_(i,j) is the jth pixel value of the ith image (each row is an image).
    B is a 1xM matrix, s.t. F_(1,j) is the jth pixel value of the mean image.
    S is the covariance matrix (NxN), computed as A x A^T.
    V is a NxN matrix, that contains eigenvectors of S as rows.
    U is a NxM matrix, that contains eigenfaces as rows.
    W is a NxN matrix, s.t W_(i,j) is the weight of the ith eigenface in the jth image.
    X is a 1xM matrix of the image to recognize.
    Wx is a Nx1 matrix, s.t. W_(i,0) is the weight of the ith eigenface in the X.
    */
    Matrix A = Matrix(N, M, read_training_data());
    Matrix im = A.getRow(N - 2);
    write_pgm("sample.pgm", &im);
    Matrix B = Matrix(1, M);
    /* Find the mean image */
    for (int c = 0; c < M; ++c)
    {
        double sum = 0;
        for (int r = 0; r < N; ++r)
        {
            sum += A.array[r][c];
        }
        B.array[0][c] = sum/N;
    }
    /* Output the mean image */
    write_pgm("meanimage.pgm", &B);
    /* Subtract the mean from each image */
    for (int r = 0; r < N; ++r)
    {
        for (int c = 0; c < M; ++c)
        {
            A.array[r][c] -= B.array[0][c];
            if (A.array[r][c] < 0)
            {
                A.array[r][c] = 0;
            }
        }
    }
    /* Output the normalized images */
    for (int i = 0; i < N; ++i)
    {
        Matrix image = A.getRow(i);
        std::ostringstream filename;
        filename << "normalized" << i << ".pgm";
        write_pgm(filename.str(), &image);
    }
    /* Find the covariance matrix */
    Matrix S = A*A.transpose();
    /* Find eigenvectors of the covariance matrix */
    Matrix V = eigensystem(&S).second.transpose();
    /* Find eigenfaces */
    Matrix U = Matrix(Eigenfaces, M);
    for (int r = 0; r < Eigenfaces; ++r)
    {
        Matrix eigenface = V.getRow(r)*A;

        U.array[r] = eigenface.array[0];
        double norm = 0;
        for (int i = 0; i < U.columns; i++) {
            norm += pow(U.array[r][i], 2);
        }
        norm = sqrt(norm);
        for (int i = 0; i < U.columns; i++) {
            U.array[r][i] /= norm;
        }
        /* Output eigenface */
        eigenface = scale(U.getRow(r));
        std::ostringstream filename;
        filename << "eigenface" << r << ".pgm";
        write_pgm(filename.str(), &eigenface);
    }
    /* Find weights */
    Matrix W = Matrix(Eigenfaces, N);
    for (int r = 0; r < Eigenfaces; ++r)
    {
        for (int c = 0; c < N; ++c)
        {
            W.array[r][c] = (U.getRow(r)*A.getRow(c).transpose()).array[0][0];
        }
    }
    /* Perform recognition */
    double accuracy = 0;
    for (int i = 1; i <= N; ++i)
    {
        std::stringstream filename;
        filename << Data_path << "s" << i << "/" << 10 << ".pgm";
        std::ifstream image(filename.str().c_str());
        std::vector< std::vector<double> > array;

        if (image.is_open()) {
            array.push_back(read_pgm(image));
            image.close();
        } else {
            std::cout << "Image was not opened.";
        }
        Matrix X = Matrix(1, M, array);
        /* Subtract the mean image */
        for (int c = 0; c < M; ++c)
        {
            X.array[0][c] -= B.array[0][c];
            if (X.array[0][c] < 0)
            {
                X.array[0][c] = 0;
            }
        }
        /* Find weights */
        Matrix Wx = Matrix(Eigenfaces, 1);
        for (int r = 0; r < Eigenfaces; ++r)
        {
            Wx.array[r][0] = (U.getRow(r)*X.transpose()).array[0][0];
        }
        /* Find the closest face from the trainig set */
        double min_distance = 0;
        int image_number = 0;
        for (int image = 0; image < N; ++image)
        {
            double distance = 0;
            for (int eigenface = 0; eigenface < Eigenfaces; ++eigenface)
            {
                distance += fabs(W.array[eigenface][image] - Wx.array[eigenface][0]);
            }
            if (distance < min_distance || image == 0)
            {
                min_distance = distance;
                image_number = image;
            }
        }
        std::cout << i << ". " << image_number + 1 << " " << min_distance << std::endl;
        if (i == image_number + 1)
        {
            accuracy = accuracy + 1;
        }
    }
    std::cout << accuracy/N;
    return 0;
}

