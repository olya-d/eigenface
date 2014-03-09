#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

using namespace std;

const int Faces = 2;
const int Samples = 1;
const int Width = 92;
const int Height = 112;
const string Data_path = "/Users/olga_andreyeva/Dropbox/Developer/eigenface/faces/";
const int N = Faces*Samples; // Total number of images
const int M = Width*Height; // Resolution


struct Matrix {
    int rows;
    int columns;
    double **array;
};


void print_matrix(Matrix* m) {
    cout << "[" << endl;
    for (int r = 0; r < m->rows; ++r)
    {
        cout << "[";
        for (int c = 0; c < m->columns; ++c)
        {
            cout << m->array[r][c];
            if (c != m->columns - 1)
            {
                cout << ", ";
            }
        }
        cout << "]";
        if (r != m->rows - 1)
        {
            cout << "," << endl;
        }
    }
    cout << endl << "]" << endl;
}


double* read_pgm(ifstream& file, int size=M) {
    double* values = new double[size];
    string line;
    getline(file, line); // Skip P2 line
    getline(file, line); // Skip width line
    getline(file, line); // Skip height line
    getline(file, line); // Skip max value line

    int val;
    int count = 0;
    while(file >> val)
    {
        values[count] = val;
        ++count;
    }
    return values;
}

double** read_training_data() {
    /*
    Returns pointer to the NxM array a, s.t
    a[i][j] is jth value of ith image.
    */
    double **array = new double*[N];
    for (int face = 0; face < Faces; ++face)
    {
        array[face] = new double[M];
        for (int sample = 0; sample < Samples; ++sample)
        {
            ostringstream filename;
            filename << Data_path << "s" << face + 1 << "/" << sample  + 1 << ".pgm";
            ifstream image(filename.str().c_str());

            if (image.is_open()) {
                array[face*Samples + sample] = read_pgm(image);
                image.close();
            } else {
                cout << "Image was not opened.";
            }
        }
    }
    return array;
}


Matrix transpose(Matrix *m) {
    Matrix result;
    result.rows = m->columns;
    result.columns = m->rows;
    result.array = new double*[result.rows];
    for (int r = 0; r < result.rows; ++r)
    {
        result.array[r] = new double[result.columns];
        for (int c = 0; c < result.columns; ++c)
        {
            result.array[r][c] = m->array[c][r];
        }
    }
    return result;
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

    return 0;
}

