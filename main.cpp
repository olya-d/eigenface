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


struct Matrix {
    int rows;
    int columns;
    std::vector< std::vector<double> > array;
};


void print_matrix(Matrix* m) {
    std::cout << "[" << std::endl;
    for (int r = 0; r < m->rows; ++r)
    {
        std::cout << "[";
        for (int c = 0; c < m->columns; ++c)
        {
            std::cout << m->array[r][c];
            if (c != m->columns - 1)
            {
                std::cout << ", ";
            }
        }
        std::cout << "]";
        if (r != m->rows - 1)
        {
            std::cout << "," << std::endl;
        }
    }
    std::cout << std::endl << "]" << std::endl;
}


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


int main(int argc, const char * argv[])
{
    // First create matrix with images as rows
    Matrix images;
    images.rows = N;
    images.columns = M;
    images.array = read_training_data();
    // Then transpose (images as columns)
    images = transpose(&images);
    // Normalize images by subtracting the mean
    Matrix mean_image = mean_column(&images);
    subtract_from_columns(&images, &mean_image);

    return 0;
}

