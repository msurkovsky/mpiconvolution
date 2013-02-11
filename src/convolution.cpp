#include "mpi.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define TAG_BLOCK        0
#define TAG_BLOCK_WIDTH  1
#define TAG_BLOCK_HEIGHT 2

float *generate_matrix(int, int, float);
float *convolve(float *, int, int, float *, int, int);
void print_matrix(float*, int, int);

int main(int argc, char *argv[]) {

    int numtasks, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    MPI_Status stats[numtasks-1];
    MPI_Request reqs[numtasks-1];

    if (argc < 3) {
        printf("You must enter dividing of a grid!\n");
        return 0;
    }

    int mask_w = 5, mask_h = 5;
    float mask[] = {
        0.0, 1.0, 0.0,
        1.0, 1.0, 1.0,
        0.0, 1.0, 0.0
    };
    int mw_overlay = mask_w / 2;
    int mh_overlay = mask_h / 2;

    int matrix_w = 24, matrix_h = 24;

    int x_size = atoi(argv[1]);
    int y_size = atoi(argv[2]);

    int original_block_width = matrix_w / x_size;
    int original_block_height = matrix_h / y_size;
    int block_width = original_block_width + 2 * mw_overlay;
    int block_height = original_block_height + 2 * mh_overlay;
    int block_size = block_width * block_height;

    float block[block_size];
//    memset(block, 0, sizeof(int) * block_size); // set a block to the zero values.

    if (rank == 0) { // divide matrix
        float *matrix = generate_matrix(matrix_w, matrix_h, 1.0);

        if (x_size * y_size != numtasks) {
            printf("The number of processes does not correspond with dividing of grid!\n");
            MPI_Finalize();
            return 0;
        }

        if (matrix_w / x_size < mask_w || matrix_h / y_size < mask_h) {
            printf("W: %d / %d < %d\n", matrix_w, x_size, mask_w);
            printf("H: %d / %d < %d\n", matrix_h, y_size, mask_h);
            printf("This example is not suitable for parallel processing!\n");
            MPI_Finalize();
            return 0;
        }

        // divide matrix into smaller pieces
        float blocks[numtasks][block_size];
        // initialize the blocks to the zero values
        memset(blocks, 0, sizeof(float) * numtasks * block_size);
        // fill block
        for (int i = 0; i < matrix_w; i++) {
            for (int j = 0; j < matrix_h; j++) {
                int x = i / original_block_width;
                int y = j / original_block_height;
                blocks[y * x_size + x][((j % original_block_height) + mh_overlay) * block_width + 
                                        (i % original_block_width)  + mw_overlay] =
                    matrix[j * matrix_w + i];
            }
        }

        print_matrix(matrix, matrix_w, matrix_h);
        for (int n = 0; n < y_size; n++) {
            int j1 = n * original_block_height;         // from top
            int j2 = j1 + original_block_height - 1;    // from bottom
            printf("{%d, %d}\n", j1, j2);
            for (int m = 1; m <= mh_overlay; m++) {
                int mm = mh_overlay - m;
                if (j1 - m >= 0) {
                    printf("j1-m = (%d-%d) = %d\n", j1, m, j1-m);
                    // copy from left-hand side
                    int y = j1 / original_block_height; // i need to get the right block (process);
                    for (int i = 0; i < matrix_w; i++) {
                        int x = i / original_block_width;
                        printf("proc(%d):[%d, %d] = %6.2f\n", 
                                y * x_size + x,
                                i % original_block_width + mw_overlay,
                                mm,
                                matrix[(j1-m) * matrix_w + i]); 
                        blocks[y * x_size + x][mm * block_width + 
                                               (i % original_block_width) + mw_overlay] =
                            matrix[(j1-m) * matrix_w + i];
                    }
                }


                if (j2 + m < matrix_h) {
                    // copy from right-hand side
                }
            }
        }
        // zero process
//        memcpy(block, blocks[0], sizeof(int) * block_size);
//
//        // sent to others processes
//        for (int t = 1; t < numtasks; t++) {
//            MPI_Isend(&blocks[t], block_size, MPI_FLOAT, t, TAG_BLOCK, MPI_COMM_WORLD, &reqs[t-1]);
//        }
//        MPI_Waitall(numtasks-1, reqs, stats);

        for (int t = 0; t < numtasks; t++) {
            printf("proc: %i:\n", t);
            print_matrix(blocks[t], block_width, block_height);
            printf("\n");
        }
        delete [] matrix;

    } else {

//        MPI_Irecv(&block, block_size, MPI_FLOAT, 0, TAG_BLOCK, MPI_COMM_WORLD, &reqs[rank-1]);
//        MPI_Wait(&reqs[rank-1], &stats[rank-1]);
    }

//    // it works every process alone
//    float *conv = convolve(block, block_width, block_height, mask, mask_w, mask_h);
//
//    // print a specific part on each process
//    for (int i = 0; i < block_width; i++) {
//        for (int j = 0; j < block_height; j++) {
//            printf("%.2f ", conv[j * block_width + i]);
//        }
//        printf("\n");
//    }
//    printf("\n");
//
//    delete [] conv;
    MPI_Finalize();
    return 0;
}

float *generate_matrix(int width, int height, float value) {
    float *out = new float[width * height];

    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            out[height * j + i] = value++;
        }
    }

    return out;
}

float *convolve(float *grid, int g_width, int g_height, 
                float *mask, int m_width, int m_height) {

    float *out = new float[g_width * g_height];

    int mw_half = m_width / 2;
    int mh_half = m_height / 2;

    int x, y;
    for (int i = 0; i < g_width; i++) {
        for (int j = 0; j < g_height; j++) {
            float sum = 0.0;
            int count = 0;
            for (int k = -mw_half; k <= mw_half; k++) {
                for (int l = -mh_half; l <= mh_half; l++) {
                    x = i + k;
                    y = j + l;
                    if (x >= 0 && x < g_width && y >= 0 && y < g_height) {
                        sum += (grid[g_height * y + x] * 
                                mask[m_height * (mh_half+l) + (mw_half+k)]); 
                        count++;
                    }
                }
            }
            out[g_height * j + i] = sum /count;
        }
    }

    return out;
}

void print_matrix(float *matrix, int width, int height) {
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            printf("%6.2f ", matrix[j * width + i]);
        }
        printf("\n");
    }
}
