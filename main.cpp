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
const string Data_path = "/Users/olga_andreyeva/Dropbox/Developer/eigenfaces_c/faces_ascii/";
const int N = Faces*Samples; // Total number of images
const int M = Width*Height; // Resolution


struct Matrix {
    int rows;
    int columns;
    double **array;
};


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


int main(int argc, const char * argv[])
{
    Matrix images;
    images.rows = N;
    images.columns = M;
    images.array = read_training_data();

    return 0;
}

