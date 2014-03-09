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


int main(int argc, const char * argv[])
{
    Matrix images;

    return 0;
}

