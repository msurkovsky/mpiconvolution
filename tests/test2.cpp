#include "../src/matrix.h"
#include <cstdio>

int main() {

//    float m[] = {
//        -1.0, -1.0, 0.0,
//        -1.0,  0.0, 1.0,
//         0.0,  1.0, 1.0
//    };

//    float m[] =
//    {
//        -1, -1, -1, -1, -1,
//        -1,  2,  2,  2, -1,
//        -1,  2,  8,  2, -1,
//        -1,  2,  2,  2, -1,
//        -1, -1, -1, -1, -1,
//    };

//    float m[] =
//    {
//        -2, -1, 0, 1, 2,
//        -2, -1, 0, 1, 2,
//        -2, -1, 0, 1, 2,
//        -2, -1, 0, 1, 2,
//        -2, -1, 0, 1, 2,
//    };

//    float m[] = {
//        -1, 1, 0, 1, 0, 1, 0,
//        1, -1, 1, 0, 1, 0, 1,
//        0, 1, -1, 1, 0, 1, 0,
//        1, 0, 1, -1, 1, 0, 1,
//        0, 1, 0, 1, -1, 1, 0,
//        1, 0, 1, 0, 1, -1, 1,
//        0, 1, 0, 1, 0, 1, -1,
//    };

    float m[] = {
        -1, 0,  1,
         0, 1,  0,
         1, 0, -1,
    };

//    float m[] = {
//        0, -1,
//        1,  0,
//    };

    int mask_width = 3;
    int mask_height = 3;
    int x_overlay = 1;
    int y_overlay = 1;

    Matrix<float> mask(mask_width, mask_height, m);

    Matrix<float> matrix(32, 32);
    for (int i = 0; i < matrix.get_width(); i++) {
        for (int j = 0; j < matrix.get_height(); j++) {
//            matrix.set(i, j, j * matrix.get_width() + i);
            if ((i * j) % 4 == 0) {
                matrix.set(i, j, 3);
            } else {
                matrix.set(i, j, 1);
            }
        }
    }

    printf("\n");
    matrix.print(1);
    Matrix<float> *matrices = matrix.divide(2, 2, x_overlay, y_overlay);
    for (int i = 0; i < 4; i++) {
        matrices[i].convolve(mask);
    }
    printf("\n");

    Matrix<float> *original = matrix.divide(1, 1, x_overlay, y_overlay);
    printf("Original:\n");
    original[0].convolve(mask);
    Matrix<float> orig_convolved = Matrix<float>::join(32, 32, 1, 1, original, x_overlay, y_overlay);
    orig_convolved.print(2);
    printf("\n");
    delete [] original;

    Matrix<float> after = Matrix<float>::join(32, 32, 2, 2, matrices, x_overlay, y_overlay);
    printf("Divided:\n");
    after.print(2);

    delete [] matrices;
    return 0;
}
