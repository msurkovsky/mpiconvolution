#include "mpi.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define TAG_BLOCK        0
#define TAG_BLOCK_WIDTH  1
#define TAG_BLOCK_HEIGHT 2

struct Range {

    int low;
    int heigh;

    Range() {
        this->low = 0;
        this->heigh = 0;
    }

    Range(int low, int heigh) {
        this->low = low;
        this->heigh = heigh;
    }
};

float *generate_matrix(int, int, float);
float *convolve(const float *, int, int, const float *, int, int);
void print_matrix(const float*, int, int);

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

    int mask_w = 3, mask_h = 3;
    const float mask[] = {
        0.0, 1.0, 0.0,
        1.0, 1.0, 1.0,
        0.0, 1.0, 0.0
    };

    int mw_overlay = mask_w / 2;
    int mh_overlay = mask_h / 2;

    int matrix_w = 12, matrix_h = 12;

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

        print_matrix(matrix, matrix_w, matrix_h);
        printf("\n");

        // divide matrix into smaller pieces
        float blocks[numtasks][block_size];
        // initialize the blocks to the zero values
        memset(blocks, 0, sizeof(float) * numtasks * block_size);
        // fill blocks

        Range x_ranges[x_size];
        Range y_ranges[y_size];
        int b,e;
        for (int i = 0; i < x_size; i++) {
           b = i * original_block_width;
           e = b + original_block_width - 1;
           x_ranges[i] = Range(b - mw_overlay, e + mw_overlay);
        }

        for (int i = 0; i < y_size; i++) {
            b = i * original_block_height;
            e = b + original_block_height - 1;
            y_ranges[i] = Range(b - mh_overlay, e + mh_overlay);
        }

        for (int y = 0; y < y_size; y++) {
            Range y_range = y_ranges[y];
            for (int j = y_range.low, j1=0; j <= y_range.heigh; j++, j1++) {

                if (j < 0 || j >= matrix_h) continue;

                for (int x = 0; x < x_size; x++) {
                    Range x_range = x_ranges[x];
                    int task = y * x_size + x;

                    if (x_range.low < 0) {
                        // shift dest pointer
                        int end_overlay = mw_overlay;
                        if (x_range.heigh >= matrix_w) {
                            end_overlay = 2 * mw_overlay;
                        }
                        memcpy(&blocks[task][j1 * block_width + mw_overlay],
                               &matrix[j * matrix_w], // there is + 0 instead of x_range.low
                               sizeof(float) * (block_width - end_overlay));
                    } else if (x_range.heigh >= matrix_w) {
                        // cut block size
                        memcpy(&blocks[task][j1 * block_width],
                               &matrix[j * matrix_w + x_range.low],
                               sizeof(float) * (block_width - mw_overlay));
                    } else {
                        memcpy(&blocks[task][j1 * block_width],
                               &matrix[j * matrix_w + x_range.low],
                               sizeof(float) * block_width);
                    }
                }
            }
        }

        for (int t = 0; t < numtasks; t++) {
//            float *conv = convolve(
//                blocks[t], block_width, block_height,
//                mask, mask_w, mask_h);

            printf("proc: %d:\n", t);
            print_matrix(blocks[t], block_width, block_height);
//            printf("\n");
//            print_matrix(conv, block_width, block_height);
            printf("\n");
//            delete [] conv;
        }

//        // zero process
//        memcpy(block, blocks[0], sizeof(int) * block_size);
//
//        // sent to others processes
//        for (int t = 1; t < numtasks; t++) {
//            MPI_Isend(&blocks[t], block_size, MPI_FLOAT, t, TAG_BLOCK, MPI_COMM_WORLD, &reqs[t-1]);
//        }
//        MPI_Waitall(numtasks-1, reqs, stats);

//        for (int t = 0; t < numtasks; t++) {
//            printf("proc: %i:\n", t);
//            print_matrix(blocks[t], block_width, block_height);
//            printf("\n");
//        }
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
            out[height * j + i] = value;
        }
    }

    return out;
}

float *convolve(const float *grid, int g_width, int g_height,
                const float *mask, int m_width, int m_height) {

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
//            out[g_height * j + i] = sum /count;
            out[g_height * j + i] = sum;
        }
    }

    return out;
}

void print_matrix(const float *matrix, int width, int height) {
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            printf("%6.2f ", matrix[j * width + i]);
        }
        printf("\n");
    }
}
