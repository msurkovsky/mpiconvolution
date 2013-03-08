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

    if (argc < 5) {
        printf("You must enter size and dividing of a grid!\n");
        return 0;
    }

    unsigned int matrix_width  = atoi(argv[1]);
    unsigned int matrix_height = atoi(argv[2]);

    unsigned int size_x  = atoi(argv[3]);
    unsigned int size_y  = atoi(argv[4]);

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

    Matrix<float> block;
    unsigned block_size = (matrix_width/size_x + 2*x_overlay) *
        (matrix_height/size_y + 2*y_overlay);

    if (rank == 0) {
        if (size_x * size_y != numtasks) {
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

        Matrix<float> matrix = generate_matrix(matrix_width, matrix_height, 2.0);

        Matrix<float> *matrices = matrix.divide(
                size_x, size_y, x_overlay, y_overlay);

        // zero process
        block = matrices[0];

        // send to others processes
        char **buffers = new char*[numtasks-1];
        for (int t = 1; t < numtasks; t++) {
            // send data
            buffers[t-1] = matrices[t].serialize();
            unsigned int size = matrices[t].get_size_of_serialized_data();
            MPI_Isend(buffers[t-1], size, MPI_CHAR, t, TAG_BLOCK, MPI_COMM_WORLD, &reqs[t-1]);
        }
        MPI_Waitall(numtasks-1, reqs, stats);

        // deinitialize data
        for (int i = 0; i < numtasks-1; i++) {
            delete [] buffers[i];
        }
        delete [] buffers;
        delete [] matrices;
    } else {
        // receive blocks
        unsigned int size = block_size * sizeof(float) + 2*sizeof(unsigned int);
        char *buff = new char[size];
        MPI_Irecv(buff, size, MPI_CHAR, 0, TAG_BLOCK, MPI_COMM_WORLD, &reqs[rank-1]);
        MPI_Wait(&reqs[rank-1], &stats[rank-1]);
        block.deserialize(buff);
        delete [] buff;
    }

    // every process convolve their part
    block.convolve(mask);

    if (rank == 0) {
        Matrix<float> *matrices = new Matrix<float>[numtasks];

        // zero process
        matrices[0] = block;

        //receive data from other processes

        char **buffers;
        buffers = new char*[numtasks-1];
        for (int t = 1; t < numtasks; t++) {
            unsigned int size = block_size * sizeof(float) + 2*sizeof(unsigned int);
            buffers[t-1] = new char[size];
            MPI_Irecv(buffers[t-1], size, MPI_CHAR, t, TAG_BLOCK, MPI_COMM_WORLD, &reqs[t-1]);
        }
        MPI_Waitall(numtasks-1, reqs, stats);

        // deserialize received data
        for (int t = 1; t < numtasks; t++) {
            matrices[t].deserialize(buffers[t-1]);
            delete [] buffers[t-1];
        }

        Matrix<float> matrix = Matrix<float>::join(
            matrix_width, matrix_height, size_x, size_y, matrices, x_overlay, y_overlay);
//        matrix.print(2);

        delete [] buffers;
        delete [] matrices;
    } else {
        char *buff = block.serialize();
        unsigned int size = block.get_size_of_serialized_data();
        MPI_Isend(buff, size, MPI_CHAR, 0, TAG_BLOCK, MPI_COMM_WORLD, &reqs[rank-1]);
//        MPI_Wait(&reqs[rank-1], &stats[rank-1]);
    }

    MPI_Finalize();
    return 0;
}
