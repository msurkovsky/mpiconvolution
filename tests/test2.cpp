#include "../src/matrix.h"
#include <cstdio>


Matrix<float> join(const unsigned int width,
               const unsigned int height,
               const unsigned int x_divide,
               const unsigned int y_divide,
               Matrix<float> * matrices) {

    int block_width = width / x_divide;
    int block_height = height / y_divide;


    Matrix<float> out(width, height);
    for (int i = 0; i < width; i++) {
        int x = i / block_width;
        for (int j = 0; j < height; j++) {
            int y = j / block_height;
            int block_idx = y * x_divide + x;

            // FIXME: bad access to blocks. Blocks are still extended.
            // TODO: find out how to add this method into matrix.h as class method (static)
            out.set(i, j, matrices[block_idx].get(i % block_width, j % block_height));
        }
    }

    return out;
}

int main() {

    float m[] = {
        0.0, 1.0, 0.0,
        1.0, 1.0, 1.0,
        0.0, 1.0, 0.0
    };

    Matrix<float> mask(3, 3, m);

    Matrix<float> matrix(24, 24);
    for (int i = 0; i < matrix.get_width(); i++) {
        for (int j = 0; j < matrix.get_height(); j++) {
            matrix.set(i, j, i + j);
        }
    }

    printf("\n");
    matrix.print(2);
    Matrix<float> *matrices = matrix.divide(12, 4, 3, 1, 1);
    for (int i = 0; i < 12; i++) {
        printf("block: %d:\n", i);
        matrices[i].convolve(mask);
        matrices[i].print(3);
        printf("\n");
    }
    printf("\n");

    printf("\n");
    Matrix<float> after = join(24, 24, 4, 3, matrices);
    after.print(3);

    delete [] matrices;
    return 0;
}
