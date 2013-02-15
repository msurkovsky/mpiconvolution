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

    if (rank == 0) { // divide matrix
        float *matrix = generate_matrix(matrix_w, matrix_h, 1.0);

        if (x_size * y_size != numtasks) {
            printf("The number of processes does not correspond with dividing of grid!\n");
            printf("%d * %d != %d\n", x_size, y_size, numtasks);
            MPI_Finalize();
            return 0;
        }

        if (matrix_w / x_size < mask_w || matrix_h / y_size < mask_h) {
            printf("This example is not suitable for parallel processing!\n");
            printf("W: %d / %d < %d\n", matrix_w, x_size, mask_w);
            printf("H: %d / %d < %d\n", matrix_h, y_size, mask_h);
            MPI_Finalize();
            return 0;
        }

        // divide matrix into smaller pieces
        float blocks[numtasks][block_size];
        // initialize the blocks to the zero values
        memset(blocks, 0, sizeof(float) * numtasks * block_size);
        // fill block
        int task = 0;

        print_matrix(matrix, matrix_w, matrix_h);
        printf("\n");

        for (int i = 0; i < matrix_w; i++) {
            for (int j = 0; j < matrix_h; j++) {
                int x = i / original_block_width;
                int y = j / original_block_height;
                int task = y * x_size + x;

                int m = (i % original_block_width) + mw_overlay;
                int n = (j % original_block_height) + mh_overlay;

                blocks[task][n * block_width + m] = matrix[j * matrix_w + i];

                // share data
                int xx = -1, yy = -1, mm, nn;
                if (i % original_block_width < mw_overlay && i - mw_overlay >= 0) {
                    xx = (i - mw_overlay) / original_block_width;
                    mm = block_width - mw_overlay + (i % original_block_width);
                } else if (i % original_block_width >= original_block_width - mw_overlay
                        && i + mw_overlay < matrix_w) {

                    xx = (i + mw_overlay) / original_block_width;
                    mm = (i % original_block_width) - original_block_width + mw_overlay;
                }

                // set x-overlay
                if (xx > -1) {
                    task = y * x_size + xx;
                    blocks[task][n * block_width + mm] = matrix[j * matrix_w + i];
                }

                if (j % original_block_height < mh_overlay && j - mh_overlay >= 0) {

                    yy = (j-mh_overlay) / original_block_height;
                    nn = block_height - mh_overlay + (j % original_block_height);
                } else if (j % original_block_height >= original_block_height - mh_overlay
                        && j + mh_overlay < matrix_h) {

                    yy = (j+mh_overlay) / original_block_height;
                    nn = (j % original_block_height) - original_block_height + mh_overlay;
                }

                // set y-overlay
                if (yy > -1) {
                    task = yy * x_size + x;
                    blocks[task][nn * block_width + m] = matrix[j * matrix_w + i];
                }

                // set x-y overlay
                if (xx > -1 && y > -1) {
                    task = yy * x_size + xx;
                    blocks[task][nn * block_width + mm] = matrix[j * matrix_w + i];
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
