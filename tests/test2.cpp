#include "../src/matrix.h"
#include <cstdio>

int main() {

//    float m[] = {
//        -1.0, -1.0, 0.0,
//        -1.0,  0.0, 1.0,
//         0.0,  1.0, 1.0
//    };

    float m[] =
    {
        -1, -1, -1, -1, -1,
        -1,  2,  2,  2, -1,
        -1,  2,  8,  2, -1,
        -1,  2,  2,  2, -1,
        -1, -1, -1, -1, -1,
    };

//    float m[] =
//    {
//        0, 0, 0, 0, 0, 0
//        0, 0, 0, 0, 0, 0
//        0, 0, 0, 0, 0, 0
//        0, 0, 0, 0, 0, 0
//        0, 0, 0, 0, 0, 0
//        0, 0, 0, 0, 0, 0
//    };

    Matrix<float> mask(5, 5, m);

    Matrix<float> matrix(32, 32);
    for (int i = 0; i < matrix.get_width(); i++) {
        for (int j = 0; j < matrix.get_height(); j++) {
//            matrix.set(i, j, j * matrix.get_width() + i);
            matrix.set(i, j, 1);
        }
    }

    printf("\n");
    matrix.print(1);
    Matrix<float> *matrices = matrix.divide(4, 2, 2, 2, 2);
    for (int i = 0; i < 4; i++) {
        matrices[i].convolve(mask);
    }
    printf("\n");

//    Matrix<float> *original = matrix.divide(1, 1, 1, 2, 2);
//    printf("Big:\n");
//    original[0].print(4);
//    printf("\n");

//    printf("Original:\n");
//    original[0].convolve(mask);
//    Matrix<float> orig_convolved = Matrix<float>::join(32, 32, 1, 1, original, 2, 2);
//    orig_convolved.print(4);
//    printf("\n");
//    delete [] original;

    Matrix<float> after = Matrix<float>::join(32, 32, 2, 2, matrices, 2, 2);
    printf("Divided:\n");
    after.print(4);

//    delete [] matrices;
    return 0;
}
