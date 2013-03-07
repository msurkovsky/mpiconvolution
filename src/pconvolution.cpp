#include "mpi.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "matrix.h"

#define TAG_BLOCK 0

Matrix<float> generate_matrix(int width, int height, float value) {
    Matrix<float> out(width, height);
    for (int i = 0; i < width; i++ ) {
        for (int j = 0; j < height; j++) {
            out.set(i, j, value);
        }
    }
    return out;
}

int main(int argc, char *argv[]) {

    int numtasks, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    // sends to the numtasks-1 processes
    MPI_Status stats[numtasks-1];
    MPI_Request reqs[numtasks-1];

    if (argc < 3) {
        printf("You must enter dividing of a grid!\n");
        return 0;
    }

    unsigned int size_x  = atoi(argv[1]);
    unsigned int size_y  = atoi(argv[2]);

    unsigned int mask_width  = 3;
    unsigned int mask_height = 3;
    const float mask[] = {
        0.0, 1.0, 0.0,
        1.0, 0.0, 1.0,
        0.0, 1.0, 0.0,
    };
    unsigned int x_overlay = mask_width / 2;
    unsigned int y_overlay = mask_height / 2;

    unsigned int matrix_width  = 2048;
    unsigned int matrix_height = 2048;

    int blocklenghts[3] = {matrix_width * matrix_height, 1, 1};
    MPI_Datatype types[3] = {MPI_FLOAT, MPI_UNSIGNED, MPI_UNSIGNED};
    MPI_Datatype mpi_matrix_type; 
    MPI_Aint offsets[3];

    Matrix<float> block;

    if (rank == 0) {
        if (size_y * size_y != numtasks) {
            printf("The number of processes does not correspond with dividing of grid!\n");
            printf("%d * %d != %d\n", size_x, size_y, numtasks);
            MPI_Finalize();
            return 0;
        }

        if (matrix_width / size_x < mask_width || 
                matrix_height / size_y < mask_height) {
            printf("This example is not suitable for parallel processing!\n");
            printf("For width:\n %d / %d < %d\n", matrix_width, size_x, mask_width);
            printf("For height:\n %d / %d < %d\n", matrix_height, size_y, mask_height);
            MPI_Finalize();
            return 0;
        }

        Matrix<float> matrix = generate_matrix(matrix_width, matrix_height, 1.0);

        Matrix<float> *matrices = matrix.divide(
                size_x, size_y, x_overlay, y_overlay);

        // zero process
        block = matrices[0];

        // send to others processes
        for (int t = 1; t < numtasks; t++) {
            // send
        }
        MPI_Waitall(numtasks-1, reqs, stats);
        delete [] matrices;
    } else {
        // receive blocks
    }

    MPI_Finalize();
    return 0;
}
