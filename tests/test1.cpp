#include <cstdio>
#include <cstdlib>
#include <cstring>


float *generate_matrix(int width, int height, float value);

int main(int argc, char *argv[]) {

    int numtasks = 1;

    int block_w = 2048;
    int block_h = 2048;
    int block_size = block_w * block_h;

    float * m = generate_matrix(block_w, block_h, 1.0);

//    float blocks[numtasks][block_size];
    float **blocks;
    blocks = new float*[numtasks];
    for (int t = 0; t < numtasks; t++) {
        blocks[t] = new float[block_size];
    }

    for (int t = 0; t < numtasks; t++) {
        for (int i = 0; i < block_size; i++) {
            blocks[t][i] = 0;
        }
    }

    for (int t = 0; t < numtasks; t++) {
        delete [] blocks[t];
    }
    delete [] blocks;
    delete [] m;
    return 0;
}

float *generate_matrix(int width, int height, float value) {
    float *out = new float[width * height];

    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            out[j * width + i] = value;
        }
    }

    return out;
}
