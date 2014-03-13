#include "eigen.h"
#include <cmath>
#include <algorithm>


Matrix identity_matrix(int dimension)
{
    Matrix I = Matrix(dimension, dimension);
    for (int r = 0; r < dimension; ++r)
    {
        for (int c = 0; c < dimension; ++c)
        {
            if (c == r)
            {
                I.array[c][r] = 1;
            }
        }
    }
    return I;
}

std::pair<std::vector<double>, Matrix> sort_by_eigenvalues(std::pair<std::vector<double>, Matrix> system)
{
    for (int i = 0; i < system.first.size() - 1; ++i)
    {
        int index = i;
        double value = system.first[i];
        for (int j = i + 1; j < system.first.size(); ++j)
        {
            if (system.first[j] > value)
            {
                index = j;
                value = system.first[j];
            }
        }
        if (index != i)
        {
            std::swap(system.first[i], system.first[index]);
            for (int r = 0; r < system.second.rows; ++r)
            {
                std::swap(system.second.array[r][i], system.second.array[r][index]);
            }
        }
    }
    return system;
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
                if (fabs(A.array[r][c]) >= max)
                {
                    max = fabs(A.array[r][c]);
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
            return sort_by_eigenvalues(result);
        }
        // Perform rotation
        double diff = A.array[l][l] - A.array[k][k];

        double t;
        if (fabs(A.array[k][l]) < fabs(diff)*1.0e-36)
        {
            t = A.array[k][l]/diff;
        }
        else
        {
            double phi = diff/(2.0*A.array[k][l]);
            t = 1.0/(fabs(phi) + sqrt(phi*phi + 1.0));
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
    return sort_by_eigenvalues(result);
}
