#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "matrix.h"

Matrix<float> generate_matrix(int width, int height, float value) {
    Matrix<float> out(width, height);
    for (int i = 0; i < width; i++ ) {
        for (int j = 0; j < height; j++) {
            out.set(i, j, value);
        }
    }
    return out;
}

int main(int argc, char *argv[]){ 

    unsigned int mask_width  = 3;
    unsigned int mask_height = 3;
    const float m[] = {
        0.0, 1.0, 0.0,
        1.0, 0.0, 1.0,
        0.0, 1.0, 0.0,
    };

    Matrix<float> mask(mask_width, mask_height, m);
    unsigned int x_overlay = mask_width / 2;
    unsigned int y_overlay = mask_height / 2;

    unsigned int matrix_width  = 2048;
    unsigned int matrix_height = 2048;

    Matrix<float> matrix = generate_matrix(matrix_width, matrix_height, 2.0);
    Matrix<float> *matrices = matrix.divide(1, 1, x_overlay, y_overlay);

    matrices[0].convolve(mask);
    matrix = Matrix<float>::join(
        matrix_width, matrix_height, 1, 1, matrices, x_overlay, y_overlay);

//    matrix.print(2);

    delete [] matrices;
    return 0;
}
